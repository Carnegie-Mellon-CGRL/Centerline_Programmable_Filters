from paraview.simple import *
import os

centerlines_path = r"C:/Users/Admin/Documents/Simvascular_Files/0-Centerlines/0-Final_Centerlines_Feb17/HP5159.vtp"
results_vtu_path = r"C:/Users/Admin/Documents/Simvascular_Files/000-Simulation_Results/HP5159/all_results_01000.vtu"
results_vtp_path = r"C:/Users/Admin/Documents/Simvascular_Files/000-Simulation_Results/HP5159/all_results_01000.vtp"

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
# 3.5) Clean filter to connect all centerline segments
# -------------------------

clean_pf1 = Clean(Input=pf1)

# Merge points within this distance
clean_pf1.Tolerance = 0.00001

# Recommended defaults (shown for clarity)
clean_pf1.ConvertLinesToPoints = 0   # leave topology intact
clean_pf1.ConvertPolysToLines = 0
clean_pf1.ConvertStripsToPolys = 0

# If you want to guarantee point merging:
clean_pf1.PointMerging = 1           # default = On

clean_pf1.UpdatePipeline()

# Re-retrieve the vtu and vtp
XMLUnstructuredGridReader1 = vtu_src  # results.vtu reader
XMLPolyDataReader2 = results_vtp_src      # results.vtp reader 

# -------------------------
# 4) Programmable Filter 2 on the centerline.  (Connects remaining line segments into continuous geometry)
# -------------------------
pf2 = ProgrammableFilter(Input=[clean_pf1])


pf2.Script = r"""
# The objective of this script is to clean at a very low threshold,
# and then connect all points that are currently listed as endpoints
# but are within a certain threshold of each other

import vtk
import numpy as np
from vtkmodules.util import numpy_support

# Input
input_data = self.GetInputDataObject(0,0) # acquires the data from the 
poly = vtk.vtkPolyData.SafeDownCast(input_data) # This should send poly data into a more pliable form

output = self.GetOutput()
output_data = self.GetOutputDataObject(0)

# Deep copy input to output to preserve geometry and existing arrays
output.DeepCopy(input_data)

try:
    locator = vtk.vtkStaticPointLocator()
except:
    locator = vtk.vtkPointLocator()
locator.SetDataSet(poly)
locator.BuildLocator()

# Geometry: reuse points; copy lines and append new ones
points = vtk.vtkPoints()
points.ShallowCopy(poly.GetPoints())

lines = vtk.vtkCellArray()
lines.DeepCopy(poly.GetLines())

num_cells_before = poly.GetNumberOfCells()

# -----Helper Functions-----
# Convenience to fetch GroupIds per-point (from incident cell)
group_arr_in = poly.GetCellData().GetArray("GroupIds")
def group_from_point(pid, default=0):
    if group_arr_in is None:
        return int(default)
    ids = vtk.vtkIdList()
    poly.GetPointCells(int(pid), ids)
    if ids.GetNumberOfIds() > 0:
        cid = ids.GetId(0)
        # Handle scalar (1-component) arrays robustly
        try:
            return int(group_arr_in.GetTuple1(cid))
        except AttributeError:
            return int(group_arr_in.GetValue(cid))
    return int(default)


def append_scalar(arr, val):
    # Append a single scalar to a 1-component array robustly.
    
    try:
        arr.InsertNextValue(val)
    except AttributeError:
        try:
            arr.InsertNextTypedTuple([val])
        except AttributeError:
            old_n = arr.GetNumberOfTuples()
            arr.SetNumberOfTuples(old_n + 1)
            try:
                arr.SetTuple1(old_n, val)
            except AttributeError:
                arr.SetTuple(old_n, [val])


#------ Step 1, Find current inlet and outlets -----#

# Get 'CenterlineIds' from CellData
id_array = poly.GetCellData().GetArray("CenterlineIds")
ids = numpy_support.vtk_to_numpy(id_array)

# Get 'GroupIds' from CellData
group_array = poly.GetCellData().GetArray("GroupIds")
group_ids = numpy_support.vtk_to_numpy(group_array)



# Build connectivity map from cells.  It is a dictionary that matches what other point every point connects to.
connectivity = {i: [] for i in range(poly.GetNumberOfPoints())}
for cell_id in range(poly.GetNumberOfCells()):
    cell = poly.GetCell(cell_id)
    num_pts = cell.GetNumberOfPoints()
    for i in range(num_pts - 1):
        p1 = cell.GetPointId(i)
        p2 = cell.GetPointId(i + 1)
        connectivity[p1].append(p2)
        connectivity[p2].append(p1)
for entries in connectivity:
    connectivity[entries] = list(set(connectivity[entries]))

# Map cell IDs to point IDs
point_ids = np.full(poly.GetNumberOfPoints(), -1)
for cell_id in range(poly.GetNumberOfCells()):
    cell = poly.GetCell(cell_id)
    for i in range(cell.GetNumberOfPoints()):
        pid = cell.GetPointId(i)
        point_ids[pid] = ids[cell_id]


# Find inlet point: first point from a cell with GroupId == 0
inlet_index = None
for cell_id in range(poly.GetNumberOfCells()):
    if group_ids[cell_id] == 0:
        cell = poly.GetCell(cell_id)
        inlet_index = cell.GetPointId(0)
        break

# Identify outlet points (one neighbor, not inlet)
outlet_indices = [i for i in range(poly.GetNumberOfPoints()) if len(connectivity[i]) <= 1 and i != inlet_index]

# Mark inlet and outlet points
inlet_array = np.zeros(poly.GetNumberOfPoints())
outlet_array = np.zeros(poly.GetNumberOfPoints())

inlet_array[inlet_index] = 1
for oid in outlet_indices:
    outlet_array[oid] = 1

# Step 2, Prepare/extend CenterlineIds 
existing_ids = poly.GetCellData().GetArray("CenterlineIds")

if existing_ids is not None and existing_ids.GetNumberOfComponents() == 1:
    new_centerline_ids = existing_ids.NewInstance()
    new_centerline_ids.DeepCopy(existing_ids)
    new_centerline_ids.SetName("CenterlineIds")
    rmin, rmax = existing_ids.GetRange()
    next_id = int(rmax) + 1
else:
    new_centerline_ids = vtk.vtkIntArray()
    new_centerline_ids.SetName("CenterlineIds")
    new_centerline_ids.SetNumberOfComponents(1)
    for i in range(num_cells_before):
        append_scalar(new_centerline_ids, int(next_id))
    next_id = int(num_cells_before)
# We will record per-new-cell values for a few known arrays
new_groupids_vals = []  # per new line cell
new_blanking_vals = []  # per new line cell (default 0)

# ----- Step 3, find each outlet and draw a sphere around it of the given radius r.  If you hit another point marked as an outlet, draw another line to it

outlet_mask = outlet_array > 0.5  # robust to float/VTK conversion
outlet_ids = np.flatnonzero(outlet_mask)

# Fetch points as an (N,3) NumPy array (fast vector access)
pts_np = numpy_support.vtk_to_numpy(poly.GetPoints().GetData())


r = 0.01  # threshold, in the same units as your data

idlist = vtk.vtkIdList()
pairs = []                         
# Optional: if you want to avoid duplicates deterministically
# seen = set()

for idx, i in enumerate(outlet_ids):
    # (a) get the coordinate
    pi = pts_np[i]  # shape (3,)

    # (b) find neighbors within the radius
    idlist.Reset()
    locator.FindPointsWithinRadius(r, pi, idlist)


    # (c) keep only *other* outlets and avoid double-counting
    for k in range(idlist.GetNumberOfIds()):
        j = idlist.GetId(k)
        if j <= i:                   # ensures each pair appears once
            continue
        if not outlet_mask[j]:       # only consider outlet-to-outlet
            continue
                
        pairs.append((i, j))

added_segments = 0

# Responsible for connecting the given points in a pair by a line
for (idA, idB) in pairs:
    if idA == idB:
        continue  # safety: skip zero-length

    # (Optional) Skip if a line already exists between them
    if idB in connectivity[idA]: continue

    line = vtk.vtkLine()
    line.GetPointIds().SetId(0, idA)
    line.GetPointIds().SetId(1, idB)
    lines.InsertNextCell(line)

    # Extend CenterlineIds (you already set up next_id earlier)
    append_scalar(new_centerline_ids, int(next_id))
    next_id += 1

    # Choose GroupIds/Blanking for this new cell (you already have helpers)
    new_groupids_vals.append(group_from_point(idA, default=0))
    new_blanking_vals.append(0)

    added_segments += 1


# === Produce output polydata ===
out_poly = vtk.vtkPolyData()
out_poly.SetPoints(points)
out_poly.SetLines(lines)

# ---- Preserve *Point Data* arrays (safe: we did not change points) ----
out_point_data = out_poly.GetPointData()
out_point_data.ShallowCopy(poly.GetPointData())

# ---- Build *Cell Data* arrays, extending all arrays to match new cell count ----
out_cell_data = out_poly.GetCellData()

num_cells_after = lines.GetNumberOfCells()
assert num_cells_after == num_cells_before + added_segments

# Helper to append tuples to a data array

def append_tuple(arr, values):
 
    #Append one tuple (or scalar) to a VTK data array in a version-agnostic way.
    #Works for vtkDataArray (numeric) and handles any #components.
    
    ncomp = arr.GetNumberOfComponents()

    # Normalize values to a sequence of length ncomp
    if ncomp == 1:
        if isinstance(values, (list, tuple)):
            val0 = values[0]
        else:
            val0 = values
        # Try the most specific methods first
        try:
            # Some typed arrays have InsertNextValue
            arr.InsertNextValue(val0)
            return
        except AttributeError:
            pass
        try:
            # Some wrappers provide InsertNextTypedTuple even for 1 comp
            arr.InsertNextTypedTuple([val0])
            return
        except AttributeError:
            pass
        # Fallback: grow and SetTuple1
        old_n = arr.GetNumberOfTuples()
        arr.SetNumberOfTuples(old_n + 1)
        try:
            arr.SetTuple1(old_n, val0)  # exists on vtkDataArray
        except AttributeError:
            # Very old fallback: SetTuple with a 1-length list
            arr.SetTuple(old_n, [val0])
        return
    else:
        # Multi-component
        if not isinstance(values, (list, tuple)) or len(values) != ncomp:
            values = [0] * ncomp

        try:
            arr.InsertNextTypedTuple(values)
            return
        except AttributeError:
            pass
        try:
            # Older VTK
            arr.InsertNextTuple(values)
            return
        except AttributeError:
            pass
        # Last resort: grow and SetTuple
        old_n = arr.GetNumberOfTuples()
        arr.SetNumberOfTuples(old_n + 1)
        arr.SetTuple(old_n, values)

# Extend every existing CellData array (except CenterlineIds, handled separately)
cd_in = poly.GetCellData()
for i in range(cd_in.GetNumberOfArrays()):
    arr_in = cd_in.GetArray(i)
    if arr_in is None:
        continue
    name = arr_in.GetName() or ""
    if name == "CenterlineIds":
        continue  # we will add our rebuilt array

    arr_out = arr_in.NewInstance()
    arr_out.DeepCopy(arr_in)
    arr_out.SetName(name)

    # Decide defaults per array
    if arr_out.IsA("vtkStringArray"):
        defaults = ["" for _ in range(added_segments)]
        # Append one-by-one
        for _ in range(added_segments):
            arr_out.InsertNextValue("")
    else:
        # Choose per-cell values:
        if name == "GroupIds":
            for val in new_groupids_vals:
                append_tuple(arr_out, int(val))
        elif name == "Blanking":
            for val in new_blanking_vals:
                append_tuple(arr_out, int(val))
        else:
            # Generic zero-fill with correct component count
            ncomp = arr_out.GetNumberOfComponents()
            zero = 0 if ncomp == 1 else [0]*ncomp
            for _ in range(added_segments):
                append_tuple(arr_out, zero)

    out_cell_data.AddArray(arr_out)

# Finally add the updated CenterlineIds
out_cell_data.AddArray(new_centerline_ids)

# Preserve active cell arrays if they existed
if cd_in.GetScalars():
    out_cell_data.SetActiveScalars(cd_in.GetScalars().GetName())
if cd_in.GetVectors():
    out_cell_data.SetActiveVectors(cd_in.GetVectors().GetName())
if cd_in.GetNormals():
    out_cell_data.SetActiveNormals(cd_in.GetNormals().GetName())
if cd_in.GetTCoords():
    out_cell_data.SetActiveTCoords(cd_in.GetTCoords().GetName())

# Finalize
output_data.ShallowCopy(out_poly)
"""

pf2.RequestInformationScript = ''
pf2.RequestUpdateExtentScript = ''
pf2.PythonPath = ''
pf2.UpdatePipeline()
# -------------------------
# 5) Programmable Filter 3 on the centerline.  (Determines inlet, outlets, and length)
# -------------------------
pf3 = ProgrammableFilter(Input=[pf2])


pf3.Script = r"""


import vtk
from vtkmodules.util import numpy_support
import numpy as np
from collections import deque

# Get input and output
input = self.GetInputDataObject(0, 0)
output = self.GetOutput()

# Deep copy input to output to preserve geometry and existing arrays
output.DeepCopy(input)

# Get 'CenterlineIds' from CellData
id_array = input.GetCellData().GetArray("CenterlineIds")
if id_array is None:
    raise RuntimeError("CenterlineIds array not found in CellData.")

ids = numpy_support.vtk_to_numpy(id_array)

# Get 'GroupIds' from CellData
group_array = input.GetCellData().GetArray("GroupIds")
if group_array is None:
    raise RuntimeError("GroupIds array not found in CellData.")

group_ids = numpy_support.vtk_to_numpy(group_array)

# Build connectivity map from cells
connectivity = {i: [] for i in range(input.GetNumberOfPoints())}
for cell_id in range(input.GetNumberOfCells()):
    cell = input.GetCell(cell_id)
    num_pts = cell.GetNumberOfPoints()
    for i in range(num_pts - 1):
        p1 = cell.GetPointId(i)
        p2 = cell.GetPointId(i + 1)
        connectivity[p1].append(p2)
        connectivity[p2].append(p1)

# Map cell IDs to point IDs
point_ids = np.full(input.GetNumberOfPoints(), -1)
for cell_id in range(input.GetNumberOfCells()):
    cell = input.GetCell(cell_id)
    for i in range(cell.GetNumberOfPoints()):
        pid = cell.GetPointId(i)
        point_ids[pid] = ids[cell_id]

# Find inlet point: first point from a cell with GroupId == 0
inlet_index = None
for cell_id in range(input.GetNumberOfCells()):
    if group_ids[cell_id] == 0:
        cell = input.GetCell(cell_id)
        inlet_index = cell.GetPointId(0)
        break

if inlet_index is None:
    raise RuntimeError("No inlet point found with GroupId == 0.")

# Identify outlet points (one neighbor, not inlet)
outlet_indices = [i for i in range(input.GetNumberOfPoints()) if len(connectivity[i]) == 1 and i != inlet_index]

# Mark inlet and outlet points
inlet_array = np.zeros(input.GetNumberOfPoints())
outlet_array = np.zeros(input.GetNumberOfPoints())

inlet_array[inlet_index] = 1
for oid in outlet_indices:
    outlet_array[oid] = 1

inlet_vtk = numpy_support.numpy_to_vtk(inlet_array)
inlet_vtk.SetName("IsInlet")
output.GetPointData().AddArray(inlet_vtk)

outlet_vtk = numpy_support.numpy_to_vtk(outlet_array)
outlet_vtk.SetName("IsOutlet")
output.GetPointData().AddArray(outlet_vtk)

# New BFS: compute distance from inlet to all reachable points
def bfs_all_distances(start, connectivity, input):
    visited = set()
    queue = deque([(start, 0.0)])  # (point_id, cumulative_distance)
    distances = np.full(input.GetNumberOfPoints(), -1.0)
    distances[start] = 0.0

    while queue:
        current, dist = queue.popleft()
        visited.add(current)

        for neighbor in connectivity[current]:
            if neighbor not in visited:
                p1 = np.array(input.GetPoint(current))
                p2 = np.array(input.GetPoint(neighbor))
                edge_length = np.linalg.norm(p2 - p1)
                new_dist = dist + edge_length

                if distances[neighbor] < 0 or new_dist < distances[neighbor]:
                    distances[neighbor] = new_dist
                    queue.append((neighbor, new_dist))

    return distances

# Compute segment distances from inlet
segment_distance = bfs_all_distances(inlet_index, connectivity, input)

# Compute straight-line distance from inlet
straight_line_distance = np.zeros(input.GetNumberOfPoints())
inlet_coord = np.array(input.GetPoint(inlet_index))

for i in range(input.GetNumberOfPoints()):
    pt = np.array(input.GetPoint(i))
    straight_line_distance[i] = np.linalg.norm(pt - inlet_coord)

linear_vtk = numpy_support.numpy_to_vtk(straight_line_distance)
linear_vtk.SetName("LinearDistanceFromInlet")
output.GetPointData().AddArray(linear_vtk)

# Compute tortuosity only for reachable points
tortuosity = np.full(input.GetNumberOfPoints(), -1.0)
for i in range(input.GetNumberOfPoints()):
    seg_dist = segment_distance[i]
    straight_dist = straight_line_distance[i]
    if seg_dist >= 0 and straight_dist > 0:
        tortuosity[i] = seg_dist / straight_dist
tortuosity_vtk = numpy_support.numpy_to_vtk(tortuosity)
tortuosity_vtk.SetName("Tortuosity")
output.GetPointData().AddArray(tortuosity_vtk)

segment_vtk = numpy_support.numpy_to_vtk(segment_distance)
segment_vtk.SetName("CenterlineLength")
output.GetPointData().AddArray(segment_vtk)

# Summary
print(f"Inlet point index (GroupId == 0): {inlet_index}")
print(f"Number of outlet points: {len(outlet_indices)}")
print("CenterlineLength and Tortuosity arrays added.")


"""

pf3.RequestInformationScript = ''
pf3.RequestUpdateExtentScript = ''
pf3.PythonPath = ''
pf3.UpdatePipeline()

# -------------------------
# 6) Programmable Filter 4 on the results of programmable filter 3
# This filter identifies bifurcation points in the centerlines
# -------------------------
pf4 = ProgrammableFilter(Input=[pf3])


pf4.Script = r"""

import vtk
from vtkmodules.util import numpy_support
import numpy as np
from collections import deque

# Get input and output
input = self.GetInputDataObject(0, 0)
output = self.GetOutput()
output.DeepCopy(input)

# Get required arrays
id_array = input.GetCellData().GetArray("CenterlineIds")
group_array = input.GetCellData().GetArray("GroupIds")
blanking_array = input.GetCellData().GetArray("Blanking")
is_inlet_array = input.GetPointData().GetArray("IsInlet")
is_outlet_array = input.GetPointData().GetArray("IsOutlet")

if id_array is None or group_array is None or blanking_array is None or is_inlet_array is None or is_outlet_array is None:
    raise RuntimeError("Required arrays not found in CellData or PointData.")

ids = numpy_support.vtk_to_numpy(id_array)
group_ids = numpy_support.vtk_to_numpy(group_array)
blanking = numpy_support.vtk_to_numpy(blanking_array)
is_inlet = numpy_support.vtk_to_numpy(is_inlet_array)
is_outlet = numpy_support.vtk_to_numpy(is_outlet_array)

num_points = input.GetNumberOfPoints()

# Convert blanking to PointData
point_blanking = np.zeros(num_points)
point_counts = np.zeros(num_points)

for cell_id in range(input.GetNumberOfCells()):
    cell = input.GetCell(cell_id)
    for i in range(cell.GetNumberOfPoints()):
        pid = cell.GetPointId(i)
        point_blanking[pid] += blanking[cell_id]
        point_counts[pid] += 1

point_counts[point_counts == 0] = 1
point_blanking = point_blanking / point_counts
point_blanking = (point_blanking > 0.5).astype(int)

# Build connectivity map
connectivity = {i: [] for i in range(num_points)}
for cell_id in range(input.GetNumberOfCells()):
    cell = input.GetCell(cell_id)
    for i in range(cell.GetNumberOfPoints() - 1):
        p1 = cell.GetPointId(i)
        p2 = cell.GetPointId(i + 1)
        connectivity[p1].append(p2)
        connectivity[p2].append(p1)

# Compute NumConnections array
num_connections = np.array([len(connectivity[i]) for i in range(num_points)])

# Find inlet point (GroupId == 0)
inlet_index = None
for cell_id in range(input.GetNumberOfCells()):
    if group_ids[cell_id] == 0:
        inlet_index = input.GetCell(cell_id).GetPointId(0)
        break
if inlet_index is None:
    raise RuntimeError("No inlet point found.")

# Identify outlet points
outlet_indices = [i for i in range(num_points) if len(connectivity[i]) == 1 and i != inlet_index]

# BFS path finder
def bfs_path(start, end, connectivity):
    visited = set()
    queue = deque([[start]])
    while queue:
        path = queue.popleft()
        node = path[-1]
        if node == end:
            return path
        if node not in visited:
            visited.add(node)
            for neighbor in connectivity[node]:
                if neighbor not in path:
                    new_path = list(path)
                    new_path.append(neighbor)
                    queue.append(new_path)
    return []

# Initialize global cap array
global_cap_array = np.zeros(num_points)

# Process each path
for outlet in outlet_indices:
    path = bfs_path(inlet_index, outlet, connectivity)

    # Cap detection based on NumConnections change
    for j in range(1, len(path)):
        current = path[j]
        previous = path[j - 1]
        if num_connections[current] != num_connections[previous]:
            global_cap_array[previous] = 1

    # Add point immediately after blanking region as cap
    for j in range(1, len(path)):
        prev = path[j - 1]
        curr = path[j]
        if point_blanking[prev] == 1 and point_blanking[curr] == 0:
            global_cap_array[curr] = 1

# Remove cap points based on updated criteria
for idx in range(num_points):
    if global_cap_array[idx] == 1:
        if num_connections[idx] <= 2 or point_blanking[idx] == 1 or is_inlet[idx] == 1:
            global_cap_array[idx] = 0

# Add final combined cap array
vtk_cap_array = numpy_support.numpy_to_vtk(global_cap_array)
vtk_cap_array.SetName("CapPoints")
output.GetPointData().AddArray(vtk_cap_array)

print(f"Found {len(outlet_indices)} outlet paths from inlet point {inlet_index}.")

"""
pf4.RequestInformationScript = ''
pf4.RequestUpdateExtentScript = ''
pf4.PythonPath = ''
pf4.UpdatePipeline()

# -------------------------
# 7) Group the result of PF4 with the original VTU
# -------------------------
grp_pf4 = GroupDatasets(Input=[pf4, XMLUnstructuredGridReader1])
grp_pf4.UpdatePipeline()



# -------------------------
# 8) Programmable Filter 5 on grouped data
# This filter will break the centerlines into sections based on the 
# radius it identifies at bifurcation points
# -------------------------
pf5 = ProgrammableFilter(Input=[grp_pf4])
pf5.Script = r"""
import vtk
from vtkmodules.util import numpy_support
import numpy as np
from collections import deque

# ----------------------------
# Config
# ----------------------------
MIN_RADIUS_EPS = 1e-6   # final clamp if everything else fails

# ----------------------------
# Helpers
# ----------------------------
def build_connectivity(centerline):
    npts = centerline.GetNumberOfPoints()
    conn = {i: [] for i in range(npts)}
    for cell_id in range(centerline.GetNumberOfCells()):
        cell = centerline.GetCell(cell_id)
        m = cell.GetNumberOfPoints()
        for j in range(m - 1):
            p1 = cell.GetPointId(j)
            p2 = cell.GetPointId(j + 1)
            conn[p1].append(p2)
            conn[p2].append(p1)
    return conn

def estimate_normal_at_point(pid, centerline_points, connectivity):
    
    # Estimate a slicing normal at a point using any neighbor direction.
    # Prefer the first neighbor; if none or degenerate, fallback to +X.
    
    p = centerline_points[pid]
    nbs = connectivity.get(pid, [])
    for nb in nbs:
        v = centerline_points[nb] - p
        n = np.linalg.norm(v)
        if np.isfinite(n) and n > 0:
            return (v / n).astype(float)
    return np.array([1.0, 0.0, 0.0], dtype=float)

def recompute_radius_by_slice(pid, centerline_points, surface, connectivity):
    
    # Re-slice the surface at this point using an estimated normal and return max edge distance.
    # Returns 0.0 if no intersection.
     
    p = centerline_points[pid].astype(float)
    nrm = estimate_normal_at_point(pid, centerline_points, connectivity)

    plane = vtk.vtkPlane()
    plane.SetOrigin(float(p[0]), float(p[1]), float(p[2]))
    plane.SetNormal(float(nrm[0]), float(nrm[1]), float(nrm[2]))

    cutter = vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputData(surface)
    cutter.Update()

    # Cleaning helps some slicers that produce duplicate points
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputData(cutter.GetOutput())
    cleaner.Update()

    region_filter = vtk.vtkPolyDataConnectivityFilter()
    region_filter.SetInputData(cleaner.GetOutput())
    region_filter.SetExtractionModeToClosestPointRegion()
    region_filter.SetClosestPoint(float(p[0]), float(p[1]), float(p[2]))
    region_filter.Update()

    poly = region_filter.GetOutput()
    pts = poly.GetPoints()
    if pts is None or pts.GetNumberOfPoints() == 0:
        return 0.0

    pts_np = numpy_support.vtk_to_numpy(pts.GetData())
    dists = np.linalg.norm(pts_np - p, axis=1) if pts_np.size > 0 else np.array([0.0])
    return float(np.max(dists)) if dists.size > 0 else 0.0

def average_neighbors_nonzero(pid, filtered_radius_array, connectivity):
    
    #Average immediate neighbors that are > 0. If none, expand outward BFS until some found.
    
    # Try 1-hop
    neigh = connectivity.get(pid, [])
    vals = [filtered_radius_array[n] for n in neigh if filtered_radius_array[n] > 0]
    if len(vals) > 0:
        return float(np.mean(vals))

    # BFS outward for first shell with non-zero entries
    visited = set([pid])
    q = deque([(pid, 0)])
    found = []
    current_depth = 0
    while q and len(found) == 0 and current_depth < 10:  # limit expansion to avoid huge searches
        node, depth = q.popleft()
        if depth > current_depth:
            current_depth = depth
        for nb in connectivity.get(node, []):
            if nb in visited:
                continue
            visited.add(nb)
            val = filtered_radius_array[nb]
            if val > 0:
                found.append(val)
            q.append((nb, depth + 1))
    if len(found) > 0:
        return float(np.mean(found))
    return 0.0

def find_path_iterative(connectivity, start, end):
    #Iterative DFS to find a path from start to end.
    stack = [(start, [start])]
    visited = set()

    while stack:
        current, path = stack.pop()
        if current == end:
            return path
        if current in visited:
            continue
        visited.add(current)
        for neighbor in connectivity[current]:
            if neighbor not in visited:
                stack.append((neighbor, path + [neighbor]))
    return []

# ----------------------------
# Get input blocks
# ----------------------------
grouped_input = inputs[0]
centerline = grouped_input.GetBlock(0)
surface    = grouped_input.GetBlock(1)

# Get centerline arrays
centerline_points = numpy_support.vtk_to_numpy(centerline.GetPoints().GetData())
cap_array_vtk     = centerline.GetPointData().GetArray("CapPoints")
is_inlet_vtk      = centerline.GetPointData().GetArray("IsInlet")
is_outlet_vtk     = centerline.GetPointData().GetArray("IsOutlet")

if cap_array_vtk is None or is_inlet_vtk is None or is_outlet_vtk is None:
    raise RuntimeError("Required arrays (CapPoints, IsInlet, IsOutlet) not found.")

cap_array  = numpy_support.vtk_to_numpy(cap_array_vtk)
is_inlet   = numpy_support.vtk_to_numpy(is_inlet_vtk)
is_outlet  = numpy_support.vtk_to_numpy(is_outlet_vtk)

num_points = centerline.GetNumberOfPoints()

# Build connectivity once (used by multiple steps)
connectivity = build_connectivity(centerline)

# ----------------------------
# Step 1: Compute radius at each CapPoint (map only)
# ----------------------------
radius_map = {}
for i in range(num_points):
    if cap_array[i] != 1:
        continue

    # Try to find a previous connected point by walking cells (as you had)
    prev_point = None
    for cell_id in range(centerline.GetNumberOfCells()):
        cell = centerline.GetCell(cell_id)
        ids = [cell.GetPointId(j) for j in range(cell.GetNumberOfPoints())]
        if i in ids:
            idx = ids.index(i)
            if idx > 0:
                prev_id = ids[idx - 1]
                prev_point = centerline_points[prev_id]
                break

    # If not found via cell order, use any neighbor (connectivity)
    if prev_point is None:
        nbs = connectivity.get(i, [])
        if len(nbs) > 0:
            prev_point = centerline_points[nbs[0]]
        else:
            continue

    normal = centerline_points[i] - prev_point
    norm = np.linalg.norm(normal)
    if norm == 0 or not np.isfinite(norm):
        continue
    normal = normal / norm

    # Slice the surface
    plane = vtk.vtkPlane()
    plane.SetOrigin(float(centerline_points[i][0]),
                    float(centerline_points[i][1]),
                    float(centerline_points[i][2]))
    plane.SetNormal(float(normal[0]), float(normal[1]), float(normal[2]))

    cutter = vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputData(surface)
    cutter.Update()

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputData(cutter.GetOutput())
    cleaner.Update()

    region_filter = vtk.vtkPolyDataConnectivityFilter()
    region_filter.SetInputData(cleaner.GetOutput())
    region_filter.SetExtractionModeToClosestPointRegion()
    region_filter.SetClosestPoint(float(centerline_points[i][0]),
                                  float(centerline_points[i][1]),
                                  float(centerline_points[i][2]))
    region_filter.Update()

    filtered = region_filter.GetOutput()
    points_vtk = filtered.GetPoints()
    if points_vtk is None or points_vtk.GetNumberOfPoints() == 0:
        continue

    points_np = numpy_support.vtk_to_numpy(points_vtk.GetData())
    distances = np.linalg.norm(points_np - centerline_points[i], axis=1)
    max_radius = float(np.max(distances)) if distances.size > 0 else 0.0

    radius_map[i] = max_radius

# ----------------------------
# Step 1.5: Compute inlet radius (start of traversal)
# ----------------------------
inlet_indices = np.where(is_inlet == 1)[0]
inlet_radius = None

if len(inlet_indices) > 0:
    inlet_id = int(inlet_indices[0])

    # Find next connected point (prefer via cell order; else any neighbor)
    next_point = None
    for cell_id in range(centerline.GetNumberOfCells()):
        cell = centerline.GetCell(cell_id)
        ids = [cell.GetPointId(j) for j in range(cell.GetNumberOfPoints())]
        if inlet_id in ids:
            idx = ids.index(inlet_id)
            if idx < len(ids) - 1:
                next_id = ids[idx + 1]
                next_point = centerline_points[next_id]
                break
    if next_point is None:
        nbs = connectivity.get(inlet_id, [])
        if len(nbs) > 0:
            next_point = centerline_points[nbs[0]]

    if next_point is not None:
        normal = next_point - centerline_points[inlet_id]
        norm = np.linalg.norm(normal)
        if norm != 0 and np.isfinite(norm):
            normal = normal / norm

            plane = vtk.vtkPlane()
            plane.SetOrigin(float(centerline_points[inlet_id][0]),
                            float(centerline_points[inlet_id][1]),
                            float(centerline_points[inlet_id][2]))
            plane.SetNormal(float(normal[0]), float(normal[1]), float(normal[2]))

            cutter = vtk.vtkCutter()
            cutter.SetCutFunction(plane)
            cutter.SetInputData(surface)
            cutter.Update()

            cleaner = vtk.vtkCleanPolyData()
            cleaner.SetInputData(cutter.GetOutput())
            cleaner.Update()

            region_filter = vtk.vtkPolyDataConnectivityFilter()
            region_filter.SetInputData(cleaner.GetOutput())
            region_filter.SetExtractionModeToClosestPointRegion()
            region_filter.SetClosestPoint(float(centerline_points[inlet_id][0]),
                                          float(centerline_points[inlet_id][1]),
                                          float(centerline_points[inlet_id][2]))
            region_filter.Update()

            filtered = region_filter.GetOutput()
            points_vtk = filtered.GetPoints()
            if points_vtk is not None and points_vtk.GetNumberOfPoints() > 0:
                points_np = numpy_support.vtk_to_numpy(points_vtk.GetData())
                distances = np.linalg.norm(points_np - centerline_points[inlet_id], axis=1)
                inlet_radius = float(np.max(distances)) if distances.size > 0 else None

# ----------------------------
# Step 2: Connectivity graph (already built above)
# ----------------------------
# connectivity = build_connectivity(centerline)  # already built

# ----------------------------
# Step 3: Traverse inlet-to-outlet paths and propagate filtered radius
# ----------------------------
def find_path_iterative(connectivity, start, end):
    #Iterative DFS to find a path from start to end.
    stack = [(start, [start])]
    visited = set()

    while stack:
        current, path = stack.pop()
        if current == end:
            return path
        if current in visited:
            continue
        visited.add(current)
        for neighbor in connectivity[current]:
            if neighbor not in visited:
                stack.append((neighbor, path + [neighbor]))
    return []

filtered_radius_array = np.zeros(num_points, dtype=float)

outlet_indices = np.where(is_outlet == 1)[0]
for inlet in inlet_indices:
    for outlet in outlet_indices:
        path = find_path_iterative(connectivity, int(inlet), int(outlet))
        if not path:
            continue  # Skip if no path found

        current_radius = inlet_radius  # Start with inlet radius

        for point_id in path:
            # Reassess radius at CapPoints
            if cap_array[point_id] == 1 and point_id in radius_map:
                new_radius = radius_map[point_id]
                if current_radius is None or new_radius <= 1.3 * current_radius:
                    current_radius = new_radius
                else:
                    # Skip updating radius if it grows too much
                    pass

            # Assign current radius if available
            if current_radius is not None and np.isfinite(current_radius) and current_radius > 0:
                filtered_radius_array[point_id] = float(current_radius)

# ----------------------------
# Step 4: Safeguard pass - fix zeros
# ----------------------------
initial_radius_array = filtered_radius_array.copy()
fix_code_array = np.zeros(num_points, dtype=np.int32)  # 0=OK, 1=re-slice, 2=neighbor avg, 3=global mean, 4=clamped

# 4a) Try re-slicing at each zero point
zero_idx = np.where(filtered_radius_array <= 0)[0]
for pid in zero_idx:
    r = recompute_radius_by_slice(pid, centerline_points, surface, connectivity)
    if np.isfinite(r) and r > 0:
        filtered_radius_array[pid] = float(r)
        fix_code_array[pid] = 1  # fixed by reassessment

# 4b) Try neighbor averaging for any still-zero points
zero_idx = np.where(filtered_radius_array <= 0)[0]
for pid in zero_idx:
    r = average_neighbors_nonzero(pid, filtered_radius_array, connectivity)
    if np.isfinite(r) and r > 0:
        filtered_radius_array[pid] = float(r)
        fix_code_array[pid] = 2  # fixed by neighbor average

# 4c) Global mean of non-zero radii
nonzero_vals = filtered_radius_array[filtered_radius_array > 0]
global_mean = float(np.mean(nonzero_vals)) if nonzero_vals.size > 0 else 0.0

zero_idx = np.where(filtered_radius_array <= 0)[0]
for pid in zero_idx:
    if global_mean > 0:
        filtered_radius_array[pid] = float(global_mean)
        fix_code_array[pid] = 3  # fixed by global mean

# 4d) Final clamp to epsilon
zero_idx = np.where(filtered_radius_array <= 0)[0]
for pid in zero_idx:
    filtered_radius_array[pid] = float(MIN_RADIUS_EPS)
    fix_code_array[pid] = 4  # clamped

# ----------------------------
# Step 5: Attach arrays to output
# ----------------------------
# Original (pre-fix) radius snapshot
filtered_initial_vtk = numpy_support.numpy_to_vtk(initial_radius_array)
filtered_initial_vtk.SetName("filtered_Radius_Initial")
centerline.GetPointData().AddArray(filtered_initial_vtk)

# Final filtered radius (post-fix)
filtered_vtk = numpy_support.numpy_to_vtk(filtered_radius_array)
filtered_vtk.SetName("Filtered_Radius")
centerline.GetPointData().AddArray(filtered_vtk)

# Final output
output.ShallowCopy(centerline)

"""
pf5.OutputDataSetType = 'vtkPolyData'
pf5.RequestInformationScript = ''
pf5.RequestUpdateExtentScript = ''
pf5.PythonPath = ''
pf5.UpdatePipeline()

# -------------------------
# 9) Transform arrays from cell data to point data
# -------------------------

cellToPoint = CellDatatoPointData(Input=pf5)
cellToPoint.UpdatePipeline()

# -------------------------
# 10) Group centerlines with the vtp to apply WSS
# -------------------------

grp_cellToPoint = GroupDatasets(Input=[cellToPoint, XMLPolyDataReader2])
grp_cellToPoint.UpdatePipeline()

# -------------------------
# 11) Programmable Filter 6 on grouped data
# applies WSS from vtp to centerlines
# -------------------------
pf6 = ProgrammableFilter(Input=[grp_cellToPoint])
pf6.Script = r"""
import vtk
from vtkmodules.util import numpy_support
import numpy as np
from collections import deque

# Get input blocks
grouped_input = inputs[0]
centerline = grouped_input.GetBlock(0)
surface = grouped_input.GetBlock(1)

# Get arrays
centerline_points = numpy_support.vtk_to_numpy(centerline.GetPoints().GetData())
filtered_radius_vtk = centerline.GetPointData().GetArray("Filtered_Radius")
blanking_vtk = centerline.GetPointData().GetArray("Blanking")
is_inlet_vtk = centerline.GetPointData().GetArray("IsInlet")
is_outlet_vtk = centerline.GetPointData().GetArray("IsOutlet")

if filtered_radius_vtk is None or blanking_vtk is None:
    raise RuntimeError("Required arrays (Filtered_Radius, Blanking) not found.")

filtered_radius = numpy_support.vtk_to_numpy(filtered_radius_vtk)
blanking = numpy_support.vtk_to_numpy(blanking_vtk)
is_inlet = numpy_support.vtk_to_numpy(is_inlet_vtk)
is_outlet = numpy_support.vtk_to_numpy(is_outlet_vtk)

num_points = centerline.GetNumberOfPoints()
filtered_wss_array = np.zeros(num_points)
clipping_flags_array = np.zeros(num_points, dtype=np.int8)  # New array to store use_clipping flags

# --- Build connectivity graph ---
connectivity = {i: [] for i in range(num_points)}
for cell_id in range(centerline.GetNumberOfCells()):
    cell = centerline.GetCell(cell_id)
    for j in range(cell.GetNumberOfPoints() - 1):
        p1 = cell.GetPointId(j)
        p2 = cell.GetPointId(j + 1)
        connectivity[p1].append(p2)
        connectivity[p2].append(p1)

# --- BFS path finder ---

def find_path(connectivity, start, end):
    visited = set()
    queue = deque([[start]])

    while queue:
        path = queue.popleft()
        node = path[-1]
        if node == end:
            return path
        if node not in visited:
            visited.add(node)
            for neighbor in connectivity[node]:
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    return []

# --- DFS path finder ---


def find_path_dfs(connectivity, start, end):
    stack = [(start, [start])]
    visited = set()

    while stack:
        node, path = stack.pop()
        if node == end:
            return path
        if node not in visited:
            visited.add(node)
            for neighbor in connectivity[node]:
                if neighbor not in visited:
                    stack.append((neighbor, path + [neighbor]))
    return []


# --- Collect all valid inlet-to-outlet paths ---
valid_path_points = set()
inlet_indices = np.where(is_inlet == 1)[0]
outlet_indices = np.where(is_outlet == 1)[0]

for inlet in inlet_indices:
    for outlet in outlet_indices:
        path = find_path_dfs(connectivity, inlet, outlet)
        valid_path_points.update(path)

# --- Prepare output blocks ---
output_blocks = vtk.vtkMultiBlockDataSet()
output_blocks.SetNumberOfBlocks(2)
output_blocks.SetBlock(0, centerline)

clipped_collection = vtk.vtkAppendPolyData()

for i in valid_path_points:
    point = centerline_points[i]
    radius = filtered_radius[i]
    if radius <= 0:
        radius = 1.0

    # Estimate normal direction
    if i > 0:
        prev_point = centerline_points[i - 1]
    else:
        prev_point = centerline_points[i] - np.array([1, 0, 0])

    normal = point - prev_point
    norm = np.linalg.norm(normal)
    if norm == 0 or np.isnan(norm):
        continue
    normal = normal / norm

    # Slice surface
    plane = vtk.vtkPlane()
    plane.SetOrigin(point)
    plane.SetNormal(normal)

    cutter = vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputData(surface)
    cutter.Update()

    # Apply region filter to isolate closest surface region
    region_filter = vtk.vtkPolyDataConnectivityFilter()
    region_filter.SetInputData(cutter.GetOutput())
    region_filter.SetExtractionModeToClosestPointRegion()
    region_filter.SetClosestPoint(point)
    region_filter.Update()
    filtered = region_filter.GetOutput()

    # Compute inlet_radius
    points_vtk = filtered.GetPoints()
    slice_radius = 0.0
    if points_vtk is not None and points_vtk.GetNumberOfPoints() > 0:
        points_np = numpy_support.vtk_to_numpy(points_vtk.GetData())
        distances = np.linalg.norm(points_np - point, axis=1)
        slice_radius = np.max(distances)

    # Decide whether to clip
    use_clipping = slice_radius > 1.4 * radius
    clipping_flags_array[i] = int(use_clipping)  # Store clipping flag
    if use_clipping:
        sphere = vtk.vtkSphere()
        sphere.SetCenter(point)
        sphere.SetRadius(radius)

        clipper = vtk.vtkClipPolyData()
        clipper.SetInputData(filtered)
        clipper.SetClipFunction(sphere)
        clipper.InsideOutOn()
        clipper.Update()
        clipped_slice = clipper.GetOutput()
    else:
        clipped_slice = filtered

    # Compute WSS average
    wss_data = clipped_slice.GetPointData().GetVectors("vWSS")
    if wss_data is None:
        filtered_wss_array[i] = 0.0
        continue

    magnitudes = []
    for j in range(clipped_slice.GetNumberOfPoints()):
        wss_vector = np.array(wss_data.GetTuple(j))
        mag = np.linalg.norm(wss_vector)
        if mag > 1e-8:
            magnitudes.append(mag)


    avg_wss = np.mean(magnitudes) if magnitudes else 0.0

    # Check if clipping was used and WSS is zero
    if use_clipping and avg_wss == 0.0:
        print(f"[Point {i}] WSS was zero after clipping — retrying with larger radius.")

        # Retry with larger sphere
        new_radius = 1.5 * radius
        sphere = vtk.vtkSphere()
        sphere.SetCenter(point)
        sphere.SetRadius(new_radius)

        clipper = vtk.vtkClipPolyData()
        clipper.SetInputData(filtered)
        clipper.SetClipFunction(sphere)
        clipper.InsideOutOn()
        clipper.Update()
        clipped_slice_retry = clipper.GetOutput()

        # Recompute WSS
        wss_data_retry = clipped_slice_retry.GetPointData().GetVectors("vWSS")
        if wss_data_retry is not None:
            magnitudes_retry = []
            for j in range(clipped_slice_retry.GetNumberOfPoints()):
                wss_vector = np.array(wss_data_retry.GetTuple(j))
                mag = np.linalg.norm(wss_vector)
                if mag > 1e-8:
                    magnitudes_retry.append(mag)
            avg_wss = np.mean(magnitudes_retry) if magnitudes_retry else 0.0
        filtered_radius[i] = new_radius
    filtered_wss_array[i] = avg_wss


# Attach arrays
wss_vtk = numpy_support.numpy_to_vtk(filtered_wss_array)
wss_vtk.SetName("Filtered_AverageWSS")
centerline.GetPointData().AddArray(wss_vtk)

clipping_flags_vtk = numpy_support.numpy_to_vtk(clipping_flags_array)
clipping_flags_vtk.SetName("Use_Clipping_Flag")
centerline.GetPointData().AddArray(clipping_flags_vtk)


# Final pass: fix any zero WSS values by averaging neighbors
for i in range(num_points):
    if filtered_wss_array[i] == 0.0:
        neighbors = connectivity.get(i, [])
        neighbor_values = [filtered_wss_array[n] for n in neighbors if filtered_wss_array[n] > 0.0]
        if neighbor_values:
            avg_neighbor_wss = np.mean(neighbor_values)
            filtered_wss_array[i] = avg_neighbor_wss
            print(f"[Point {i}] WSS was zero — replaced with average of neighbors: {avg_neighbor_wss}")
        else:
            print(f"[Point {i}] WSS was zero and no valid neighbors found — value remains 0.")

# Update the array with fixed values
wss_vtk_fixed = numpy_support.numpy_to_vtk(filtered_wss_array)
wss_vtk_fixed.SetName("Filtered_AverageWSS")
centerline.GetPointData().AddArray(wss_vtk_fixed)

filtered_radius_vtk_updated = numpy_support.numpy_to_vtk(filtered_radius)
filtered_radius_vtk_updated.SetName("Filtered_Radius")
centerline.GetPointData().AddArray(filtered_radius_vtk_updated)

# Output
output.ShallowCopy(centerline)

"""
pf6.OutputDataSetType = 'vtkPolyData'
pf6.RequestInformationScript = ''
pf6.RequestUpdateExtentScript = ''
pf6.PythonPath = ''
pf6.UpdatePipeline()

# -------------------------
# 12) Group centerlines with the vtu to apply Area
# -------------------------

grp_pf6 = GroupDatasets(Input=[pf6, XMLUnstructuredGridReader1])
grp_pf6.UpdatePipeline()

# -------------------------
# 13) Programmable Filter 7 on grouped data
# applies area from vtu to centerlines
# -------------------------
pf7 = ProgrammableFilter(Input=[grp_pf6])
pf7.Script = r"""

import vtk
from vtkmodules.util import numpy_support
import numpy as np
from collections import deque

# Get input blocks
grouped_input = inputs[0]
centerline = grouped_input.GetBlock(0)
surface = grouped_input.GetBlock(1)

# Get arrays
centerline_points = numpy_support.vtk_to_numpy(centerline.GetPoints().GetData())
filtered_radius_vtk = centerline.GetPointData().GetArray("Filtered_Radius")
blanking_vtk = centerline.GetPointData().GetArray("Blanking")
is_inlet_vtk = centerline.GetPointData().GetArray("IsInlet")
is_outlet_vtk = centerline.GetPointData().GetArray("IsOutlet")

if filtered_radius_vtk is None or blanking_vtk is None:
    raise RuntimeError("Required arrays (Filtered_Radius, Blanking) not found.")

filtered_radius = numpy_support.vtk_to_numpy(filtered_radius_vtk)
blanking = numpy_support.vtk_to_numpy(blanking_vtk)
is_inlet = numpy_support.vtk_to_numpy(is_inlet_vtk)
is_outlet = numpy_support.vtk_to_numpy(is_outlet_vtk)

num_points = centerline.GetNumberOfPoints()
estimated_area_array = np.zeros(num_points)
clipping_flags_array = np.zeros(num_points, dtype=np.int8)

# --- Build connectivity graph ---
connectivity = {i: [] for i in range(num_points)}
for cell_id in range(centerline.GetNumberOfCells()):
    cell = centerline.GetCell(cell_id)
    for j in range(cell.GetNumberOfPoints() - 1):
        p1 = cell.GetPointId(j)
        p2 = cell.GetPointId(j + 1)
        connectivity[p1].append(p2)
        connectivity[p2].append(p1)

# --- BFS path finder ---
def find_path(connectivity, start, end):
    visited = set()
    queue = deque([[start]])
    while queue:
        path = queue.popleft()
        node = path[-1]
        if node == end:
            return path
        if node not in visited:
            visited.add(node)
            for neighbor in connectivity[node]:
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    return []

# --- Collect all valid inlet-to-outlet paths ---
valid_path_points = set()
inlet_indices = np.where(is_inlet == 1)[0]
outlet_indices = np.where(is_outlet == 1)[0]

for inlet in inlet_indices:
    for outlet in outlet_indices:
        path = find_path(connectivity, inlet, outlet)
        valid_path_points.update(path)

# --- Prepare output blocks ---
output_blocks = vtk.vtkMultiBlockDataSet()
output_blocks.SetNumberOfBlocks(2)
output_blocks.SetBlock(0, centerline)

clipped_collection = vtk.vtkAppendPolyData()

for i in valid_path_points:
    point = centerline_points[i]
    radius = filtered_radius[i]
    if radius <= 0:
        radius = 1.0

    # Estimate normal direction
    if i > 0:
        prev_point = centerline_points[i - 1]
    else:
        prev_point = centerline_points[i] - np.array([1, 0, 0])

    normal = point - prev_point
    norm = np.linalg.norm(normal)
    if norm == 0 or np.isnan(norm):
        continue
    normal = normal / norm

    # Slice surface
    plane = vtk.vtkPlane()
    plane.SetOrigin(point)
    plane.SetNormal(normal)

    cutter = vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputData(surface)
    cutter.Update()

    # Apply region filter to isolate closest surface region
    region_filter = vtk.vtkPolyDataConnectivityFilter()
    region_filter.SetInputData(cutter.GetOutput())
    region_filter.SetExtractionModeToClosestPointRegion()
    region_filter.SetClosestPoint(point)
    region_filter.Update()
    filtered = region_filter.GetOutput()

    # Compute slice radius
    points_vtk = filtered.GetPoints()
    slice_radius = 0.0
    if points_vtk is not None and points_vtk.GetNumberOfPoints() > 0:
        points_np = numpy_support.vtk_to_numpy(points_vtk.GetData())
        distances = np.linalg.norm(points_np - point, axis=1)
        slice_radius = np.max(distances)

    # Decide whether to clip
    use_clipping = blanking[i] != 0 or slice_radius > 1.4 * radius
    clipping_flags_array[i] = int(use_clipping)
    if use_clipping:
        sphere = vtk.vtkSphere()
        sphere.SetCenter(point)
        sphere.SetRadius(radius)

        clipper = vtk.vtkClipPolyData()
        clipper.SetInputData(filtered)
        clipper.SetClipFunction(sphere)
        clipper.InsideOutOn()
        clipper.Update()
        clipped_slice = clipper.GetOutput()
    else:
        clipped_slice = filtered

    # Compute area
    triangle_filter = vtk.vtkTriangleFilter()
    triangle_filter.SetInputData(clipped_slice)
    triangle_filter.Update()

    mass_props = vtk.vtkMassProperties()
    mass_props.SetInputData(triangle_filter.GetOutput())
    clipped_area = mass_props.GetSurfaceArea()
    estimated_area_array[i] = clipped_area
    print(f"[Point {i}] Clipped area: {clipped_area}")


# Final pass: fix any zero area values by averaging neighbors
for i in range(num_points):
    if estimated_area_array[i] == 0.0:
        neighbors = connectivity.get(i, [])
        neighbor_values = [estimated_area_array[n] for n in neighbors if estimated_area_array[n] > 0.0]
        if neighbor_values:
            avg_neighbor_area = np.mean(neighbor_values)
            estimated_area_array[i] = avg_neighbor_area
            print(f"[Point {i}] Area was zero — replaced with average of neighbors: {avg_neighbor_area}")
        else:
            print(f"[Point {i}] Area was zero and no valid neighbors found — value remains 0.")


# Attach arrays
area_vtk_fixed = numpy_support.numpy_to_vtk(estimated_area_array)
area_vtk_fixed.SetName("Filtered_Area")
centerline.GetPointData().AddArray(area_vtk_fixed)

clipping_flags_vtk = numpy_support.numpy_to_vtk(clipping_flags_array)
clipping_flags_vtk.SetName("Use_Clipping_Flag")
centerline.GetPointData().AddArray(clipping_flags_vtk)

# Output
output.ShallowCopy(centerline)

"""
pf7.OutputDataSetType = 'vtkPolyData'
pf7.RequestInformationScript = ''
pf7.RequestUpdateExtentScript = ''
pf7.PythonPath = ''
pf7.UpdatePipeline()