
from paraview.simple import *
import os

centerlines_path = r"C:/Users/Admin/Documents/Simvascular_Files/0078_H_PULM_H/Models/0078_Centerlines.vtp"
results_vtu_path = r"C:/Users/Admin/Documents/Simvascular_Files/0078_H_PULM_H/Simulations/0078_Steady_Calc_2/0078_Steady_Calc_2-converted-results/all_results_02000.vtu"
results_vtp_path = r"C:/Users/Admin/Documents/Simvascular_Files/0078_H_PULM_H/Simulations/0078_Steady_Calc_2/0078_Steady_Calc_2-converted-results/all_results_02000.vtp"

# -------------------------
# 1) Load data sources
# -------------------------

if not os.path.exists(centerlines_path):
    raise RuntimeError(f"Missing file: {centerlines_path}")

if not os.path.exists(results_vtu_path):
    raise RuntimeError(f"Missing file: {results_vtu_path}")

centerlines_src = XMLPolyDataReader(FileName=[centerlines_path])
centerlines_src.UpdatePipeline()

vtu_src = XMLUnstructuredGridReader(FileName=[results_vtu_path])
vtu_src.UpdatePipeline()

if os.path.exists(results_vtp_path):
    results_vtp_src = XMLPolyDataReader(FileName=[results_vtp_path])
    results_vtp_src.UpdatePipeline()

# -------------------------
# 2) ResampleWithDataset(centerlines, vtu)
# -------------------------

resampled = ResampleWithDataset(SourceDataArrays=vtu_src, DestinationMesh=centerlines_src) 
resampled.UpdatePipeline()


# -------------------------
# 3) Programmable Filter 1 on the resampled centerlines
# -------------------------
pf1 = ProgrammableFilter(Input=[centerlines_src,resampled])

# If your Programmable Filter expects multiple inputs via 'inputs[0]' etc.,
# you can add more via the 'Input' property as a list, but you've moved to
# single-input usage with: self.GetInputDataObject(0, 0).
# Paste your "Filter 1" Python code below as a single string.

pf1.Script = r"""

import vtk
import numpy as np
from vtkmodules.util import numpy_support

# --- Inputs ---
input0 = self.GetInputDataObject(0, 0)  # reference geometry
input1 = self.GetInputDataObject(0, 1)  # dataset providing extra arrays
output = self.GetOutput()

# --- Copy geometry from input0 ---
output.DeepCopy(input0)

# --- Build locator on input1 for coordinate-based matching ---
locator = vtk.vtkPointLocator()
locator.SetDataSet(input1)
locator.BuildLocator()

tol = 1e-6  # tolerance for coordinate matching
# --- Map PointData arrays from input1 to output ---
n_points = input0.GetNumberOfPoints()
pointdata1 = input1.GetPointData()

for iarray in range(pointdata1.GetNumberOfArrays()):
    arr = pointdata1.GetArray(iarray)
    name = arr.GetName()
    ncomp = arr.GetNumberOfComponents()

    new_arr = vtk.vtkFloatArray()
    new_arr.SetName(name)
    new_arr.SetNumberOfComponents(ncomp)
    new_arr.SetNumberOfTuples(n_points)

    for pid in range(n_points):
        coord = input0.GetPoint(pid)
        match_id = locator.FindClosestPoint(coord)
        dist = np.linalg.norm(np.array(coord) - np.array(input1.GetPoint(match_id)))
        if dist < tol:
            new_arr.SetTuple(pid, arr.GetTuple(match_id))
        else:
            new_arr.SetTuple(pid, (np.nan,) * ncomp)

    output.GetPointData().AddArray(new_arr)

# --- Map CellData arrays from input1 to output ---
n_cells = input0.GetNumberOfCells()
celldata1 = input1.GetCellData()

for iarray in range(celldata1.GetNumberOfArrays()):
    arr = celldata1.GetArray(iarray)
    name = arr.GetName()
    ncomp = arr.GetNumberOfComponents()

    new_arr = vtk.vtkFloatArray()
    new_arr.SetName(name)
    new_arr.SetNumberOfComponents(ncomp)
    new_arr.SetNumberOfTuples(n_cells)

    # Simple index-based mapping (assuming same cell order)
    for cid in range(n_cells):
        if cid < input1.GetNumberOfCells():
            new_arr.SetTuple(cid, arr.GetTuple(cid))
        else:
            new_arr.SetTuple(cid, (np.nan,) * ncomp)

    output.GetCellData().AddArray(new_arr)

"""

pf1.RequestInformationScript = ''
pf1.RequestUpdateExtentScript = ''
pf1.PythonPath = ''
pf1.UpdatePipeline()


# -------------------------
# 4) Group the result of PF1 with the original VTU
# -------------------------
grp = GroupDatasets(Input=[pf1, vtu_src])
grp.UpdatePipeline()


