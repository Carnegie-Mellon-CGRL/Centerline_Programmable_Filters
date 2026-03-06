from paraview.simple import *
import os

centerlines_path = r"C:/Users/Admin/Documents/temp-2/HP5159.vtp"
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
# 2) Programmable Filter 1 on the resampled centerlines
# -------------------------


pf1 = ProgrammableFilter(Input=[centerlines_src])

# Resamples the points of the centerlines such that each data point is 0.01 cm apart.

pf1.Script = r"""

import vtk
from collections import deque

# Parameters
target_spacing = 0.01  # cm

# Input
input_data = self.GetInput()
points = input_data.GetPoints()
lines = input_data.GetLines()
point_data = input_data.GetPointData()
num_points = points.GetNumberOfPoints()
num_arrays = point_data.GetNumberOfArrays()

# Build connectivity graph
connectivity = [[] for _ in range(num_points)]
lines.InitTraversal()
id_list = vtk.vtkIdList()
while lines.GetNextCell(id_list):
    for i in range(id_list.GetNumberOfIds() - 1):
        a = id_list.GetId(i)
        b = id_list.GetId(i + 1)
        connectivity[a].append(b)
        connectivity[b].append(a)

# Find inlet and all outlets
is_inlet = point_data.GetArray("IsInlet")
is_outlet = point_data.GetArray("IsOutlet")
inlet_id = next(i for i in range(num_points) if is_inlet and is_inlet.GetValue(i) == 1)
outlet_ids = [i for i in range(num_points) if is_outlet and is_outlet.GetValue(i) == 1]

# Prepare output containers
new_points = vtk.vtkPoints()
new_lines = vtk.vtkCellArray()
new_point_data_arrays = []
for a in range(num_arrays):
    array = point_data.GetArray(a)
    new_array = vtk.vtkDoubleArray()
    new_array.SetName(array.GetName())
    new_array.SetNumberOfComponents(array.GetNumberOfComponents())
    new_point_data_arrays.append(new_array)

# Helper functions
def trace_path(inlet, outlet):
    prev = [-1] * num_points
    visited = [False] * num_points
    queue = deque()
    queue.append(inlet)
    visited[inlet] = True

    while queue:
        current = queue.popleft()
        if current == outlet:
            break
        for neighbor in connectivity[current]:
            if not visited[neighbor]:
                visited[neighbor] = True
                prev[neighbor] = current
                queue.append(neighbor)

    path = []
    current = outlet
    while current != -1:
        path.append(current)
        current = prev[current]
    path.reverse()
    return path

def interpolate_arrays(val0, val1, t):
    return [val0[c] + (val1[c] - val0[c]) * t for c in range(len(val0))]

def get_point_data(i):
    return [[point_data.GetArray(a).GetComponent(i, c)
             for c in range(point_data.GetArray(a).GetNumberOfComponents())]
            for a in range(num_arrays)]

def resample_path(path_ids):
    arc_lengths = [0.0]
    for i in range(1, len(path_ids)):
        p0 = points.GetPoint(path_ids[i - 1])
        p1 = points.GetPoint(path_ids[i])
        dist = ((p1[0] - p0[0])**2 + (p1[1] - p0[1])**2 + (p1[2] - p0[2])**2)**0.5
        arc_lengths.append(arc_lengths[-1] + dist)

    total_length = arc_lengths[-1]
    num_samples = int(total_length / target_spacing)
    sample_distances = [i * target_spacing for i in range(num_samples + 1)]

    segment_ids = []
    j = 0
    for d in sample_distances:
        while j < len(arc_lengths) - 1 and arc_lengths[j + 1] < d:
            j += 1
        t = (d - arc_lengths[j]) / (arc_lengths[j + 1] - arc_lengths[j])
        p0 = points.GetPoint(path_ids[j])
        p1 = points.GetPoint(path_ids[j + 1])
        pt = [p0[i] + t * (p1[i] - p0[i]) for i in range(3)]
        pid = new_points.InsertNextPoint(pt)

        # Interpolate PointData
        data0 = get_point_data(path_ids[j])
        data1 = get_point_data(path_ids[j + 1])
        interpolated_data = [interpolate_arrays(data0[a], data1[a], t) for a in range(num_arrays)]
        
        # Apply rounding rules
        for a, vals in enumerate(interpolated_data):
            name = new_point_data_arrays[a].GetName()

            if name == "capPoints":
                vals = [1.0 if v >= 0.5 else 0.0 for v in vals]

            elif name == "Blanking":
                vals = [1.0 if v != 0.0 else 0.0 for v in vals]

            new_point_data_arrays[a].InsertNextTuple(vals)


        segment_ids.append(pid)

    return segment_ids

# Process each outlet
for outlet_id in outlet_ids:
    path_ids = trace_path(inlet_id, outlet_id)
    segment_ids = resample_path(path_ids)

    for i in range(len(segment_ids) - 1):
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, segment_ids[i])
        line.GetPointIds().SetId(1, segment_ids[i + 1])
        new_lines.InsertNextCell(line)

# Output
output_polydata = vtk.vtkPolyData()
output_polydata.SetPoints(new_points)
output_polydata.SetLines(new_lines)
for array in new_point_data_arrays:
    output_polydata.GetPointData().AddArray(array)
output.ShallowCopy(output_polydata)

"""

pf1.RequestInformationScript = ''
pf1.RequestUpdateExtentScript = ''
pf1.PythonPath = ''
pf1.UpdatePipeline()

# -------------------------
# 3) Clean filter to connect all centerline segments
# -------------------------

clean_pf1 = Clean(Input=pf1)

# Merge points within this distance
clean_pf1.Tolerance = 0.0

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
# 4) Programmable Filter 2 on the  centerlines.
# Identifies generations of vessels and marks accordingly
# -------------------------
pf2 = ProgrammableFilter(Input=[clean_pf1])


pf2.Script = r"""

from vtkmodules.util import numpy_support
import numpy as np
from collections import deque

# Get input
input = self.GetInput()

# Get arrays
points = input.GetPoints()
num_points = input.GetNumberOfPoints()

is_inlet_vtk = input.GetPointData().GetArray("IsInlet")
if is_inlet_vtk is None:
    raise RuntimeError("Required array 'IsInlet' not found.")

is_inlet = numpy_support.vtk_to_numpy(is_inlet_vtk)

# Build connectivity graph
connectivity = {i: [] for i in range(num_points)}
for cell_id in range(input.GetNumberOfCells()):
    cell = input.GetCell(cell_id)
    for j in range(cell.GetNumberOfPoints() - 1):
        p1 = cell.GetPointId(j)
        p2 = cell.GetPointId(j + 1)
        connectivity[p1].append(p2)
        connectivity[p2].append(p1)

# Compute number of connections per point
num_connections = np.array([len(connectivity[i]) for i in range(num_points)])

# Initialize Generation array
generation_array = np.full(num_points, -1, dtype=np.int32)  # -1 means unvisited

# BFS traversal from inlet points
inlet_indices = np.where(is_inlet == 1)[0]
queue = deque()

for inlet in inlet_indices:
    generation_array[inlet] = 0
    queue.append((inlet, 0, True))  # (point_id, generation, is_first_step)

while queue:
    current, gen, is_first_step = queue.popleft()
    for neighbor in connectivity[current]:
        if generation_array[neighbor] == -1:  # not visited
            # Only increment generation if not the first step from inlet
            if is_first_step:
                next_gen = gen
            else:
                next_gen = gen + 1 if num_connections[neighbor] != num_connections[current] else gen
            generation_array[neighbor] = next_gen
            queue.append((neighbor, next_gen, False))

# Attach Generation array
generation_vtk = numpy_support.numpy_to_vtk(generation_array)
generation_vtk.SetName("Generation")
input.GetPointData().AddArray(generation_vtk)

# Output
output.ShallowCopy(input)

"""

pf2.RequestInformationScript = ''
pf2.RequestUpdateExtentScript = ''
pf2.PythonPath = ''
pf2.UpdatePipeline()

# -------------------------
# 5) Programmable Filter 3 on the results of programmable filter 2
# Rounds up extrapolated 'blanking' values
# -------------------------
pf3 = ProgrammableFilter(Input=[pf2])


pf3.Script = r"""

from vtkmodules.util import numpy_support
import numpy as np

# Get input
input_data = self.GetInput()

# Get the Blanking array from PointData
blanking_vtk = input_data.GetPointData().GetArray("Blanking")
if blanking_vtk is None:
    raise RuntimeError("Array 'Blanking' not found in PointData.")

# Convert to NumPy
blanking_np = numpy_support.vtk_to_numpy(blanking_vtk)

# Set all non-zero values to 1
normalized_blanking = np.where(blanking_np != 0, 1, 0).astype(np.int8)

# Convert back to VTK array
normalized_vtk = numpy_support.numpy_to_vtk(normalized_blanking)
normalized_vtk.SetName("blanking")

# Add to output
input_data.GetPointData().AddArray(normalized_vtk)

# Output
output.ShallowCopy(input_data)

"""
pf3.RequestInformationScript = ''
pf3.RequestUpdateExtentScript = ''
pf3.PythonPath = ''
pf3.UpdatePipeline()

# -------------------------
# 5) Programmable Filter 4 on the results of programmable filter 3
# Re-obtains the inlet and outlet points, and measures length again
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

# Deep copy input to output to preserve geometry and existing arrays
output.DeepCopy(input)

num_points = input.GetNumberOfPoints()

# Get 'CenterlineIds' from PointData
id_array = input.GetPointData().GetArray("CenterlineIds")
if id_array is None:
    raise RuntimeError("CenterlineIds array not found in PointData.")
ids = numpy_support.vtk_to_numpy(id_array)

# Get 'GroupIds' from PointData
group_array = input.GetPointData().GetArray("GroupIds")
if group_array is None:
    raise RuntimeError("GroupIds array not found in PointData.")
group_ids = numpy_support.vtk_to_numpy(group_array)

# Build connectivity map from cells
connectivity = {i: [] for i in range(num_points)}
for cell_id in range(input.GetNumberOfCells()):
    cell = input.GetCell(cell_id)
    num_pts = cell.GetNumberOfPoints()
    for i in range(num_pts - 1):
        p1 = cell.GetPointId(i)
        p2 = cell.GetPointId(i + 1)
        connectivity[p1].append(p2)
        connectivity[p2].append(p1)

# Find inlet point: first point with GroupId == 0
inlet_index = None
for pid in range(num_points):
    if group_ids[pid] == 0:
        inlet_index = pid
        break

if inlet_index is None:
    raise RuntimeError("No inlet point found with GroupId == 0.")

# Identify outlet points (one neighbor, not inlet)
outlet_indices = [i for i in range(num_points) if len(connectivity[i]) == 1 and i != inlet_index]

# Mark inlet and outlet points
inlet_array = np.zeros(num_points)
outlet_array = np.zeros(num_points)

inlet_array[inlet_index] = 1
for oid in outlet_indices:
    outlet_array[oid] = 1

inlet_vtk = numpy_support.numpy_to_vtk(inlet_array)
inlet_vtk.SetName("IsInlet")
output.GetPointData().AddArray(inlet_vtk)

outlet_vtk = numpy_support.numpy_to_vtk(outlet_array)
outlet_vtk.SetName("IsOutlet")
output.GetPointData().AddArray(outlet_vtk)

# BFS: compute distance from inlet to all reachable points
def bfs_all_distances(start, connectivity, input):
    visited = set()
    queue = deque([(start, 0.0)])
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
straight_line_distance = np.zeros(num_points)
inlet_coord = np.array(input.GetPoint(inlet_index))

for i in range(num_points):
    pt = np.array(input.GetPoint(i))
    straight_line_distance[i] = np.linalg.norm(pt - inlet_coord)

linear_vtk = numpy_support.numpy_to_vtk(straight_line_distance)
linear_vtk.SetName("LinearDistanceFromInlet")
output.GetPointData().AddArray(linear_vtk)

# Compute tortuosity only for reachable points
tortuosity = np.full(num_points, -1.0)
for i in range(num_points):
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
pf4.RequestInformationScript = ''
pf4.RequestUpdateExtentScript = ''
pf4.PythonPath = ''
pf4.UpdatePipeline()

# -------------------------
# 6) Calculator to obtain pressure in mmHg
# -------------------------


calculator = Calculator(Input=pf4)
calculator.ResultArrayName = 'Pressure (mmHg)'
calculator.Function = 'pressure / 1333.22'
calculator.UpdatePipeline()

# -------------------------
# 7) Programmable Filter 5 
# Calculates Flow, 'Power', and Resistance by Ohm's Law
# -------------------------
pf5 = ProgrammableFilter(Input=[calculator])
pf5.Script = r"""

import vtk
from vtkmodules.util import numpy_support
import numpy as np
from collections import deque

# -------------------------------------------------------------------------
# Input / Output (as requested)
# -------------------------------------------------------------------------
input = self.GetInputDataObject(0, 0)
output = self.GetOutput()
output.DeepCopy(input)

poly = input
num_points = poly.GetNumberOfPoints()

# -------------------------------------------------------------------------
# Required Point Data Arrays
# -------------------------------------------------------------------------
pd = poly.GetPointData()

filtered_area_vtk = pd.GetArray("Filtered_Area")
blanking_vtk = pd.GetArray("Blanking")
is_inlet_vtk = pd.GetArray("IsInlet")
is_outlet_vtk = pd.GetArray("IsOutlet")
velocity_vtk = pd.GetArray("velocity")
pressure_vtk = pd.GetArray("pressure")

if (filtered_area_vtk is None or blanking_vtk is None or
    velocity_vtk is None or pressure_vtk is None or
    is_inlet_vtk is None or is_outlet_vtk is None):
    raise RuntimeError("Missing required Point Data arrays.")

# -------------------------------------------------------------------------
# Convert to NumPy
# -------------------------------------------------------------------------
points_np     = numpy_support.vtk_to_numpy(poly.GetPoints().GetData())
filtered_area = numpy_support.vtk_to_numpy(filtered_area_vtk)
blanking      = numpy_support.vtk_to_numpy(blanking_vtk)
is_inlet      = numpy_support.vtk_to_numpy(is_inlet_vtk)
is_outlet     = numpy_support.vtk_to_numpy(is_outlet_vtk)
velocity_np   = numpy_support.vtk_to_numpy(velocity_vtk)
pressure_np   = numpy_support.vtk_to_numpy(pressure_vtk)

# -------------------------------------------------------------------------
# Allocate Output Arrays
# -------------------------------------------------------------------------
flow_array       = np.zeros(num_points)
power_array      = np.zeros(num_points)
resistance_array = np.zeros(num_points)
normals_array    = np.zeros((num_points, 3))

# -------------------------------------------------------------------------
# Build Connectivity Graph
# -------------------------------------------------------------------------
connectivity = {i: [] for i in range(num_points)}

for cell_id in range(poly.GetNumberOfCells()):
    cell = poly.GetCell(cell_id)
    for j in range(cell.GetNumberOfPoints() - 1):
        p1 = cell.GetPointId(j)
        p2 = cell.GetPointId(j + 1)
        connectivity[p1].append(p2)
        connectivity[p2].append(p1)

# -------------------------------------------------------------------------
# BFS for inlet → outlet paths
# -------------------------------------------------------------------------
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
                queue.append(path + [neighbor])

    return []

# -------------------------------------------------------------------------
# Identify Valid Path Points
# -------------------------------------------------------------------------
valid_path_points = set()
inlet_indices  = np.where(is_inlet == 1)[0]
outlet_indices = np.where(is_outlet == 1)[0]

for inlet in inlet_indices:
    for outlet in outlet_indices:
        path = find_path(connectivity, inlet, outlet)
        valid_path_points.update(path)

# -------------------------------------------------------------------------
# Compute Flow, Power, Resistance
# -------------------------------------------------------------------------
for i in valid_path_points:

    # Estimate normal
    if i > 0:
        prev_point = points_np[i - 1]
    else:
        prev_point = points_np[i] - np.array([1.0, 0.0, 0.0])

    tangent = points_np[i] - prev_point
    norm = np.linalg.norm(tangent)
    if norm == 0.0 or np.isnan(norm):
        continue

    normal = tangent / norm
    normals_array[i] = normal

    # Dot velocity with normal
    v_vec = velocity_np[i]
    v_normal = np.dot(v_vec, normal)

    A = filtered_area[i]
    if A <= 0.0:
        continue

    flow = abs(v_normal * A)
    flow_array[i] = flow

    if abs(flow) > 1e-8:
        power_array[i] = pressure_np[i] * flow
        resistance_array[i] = pressure_np[i] / flow

# -------------------------------------------------------------------------
# >>> SAFETY FIX: FLOW JUMP STABILIZATION <<<
# -------------------------------------------------------------------------

FLOW_JUMP_FACTOR = 10.0  # max allowed ratio between neighbors

fixed_flow = flow_array.copy()

for i in valid_path_points:

    neighbors = connectivity[i]
    neigh_flows = [flow_array[n] for n in neighbors if n in valid_path_points]

    if len(neigh_flows) == 0:
        continue

    fi = flow_array[i]

    # If any zero neighbors, handle separately
    if fi <= 0 or min(neigh_flows) <= 0:
        nz = [f for f in neigh_flows if f > 0]
        if len(nz) > 0:
            fixed_flow[i] = sum(nz) / len(nz)
        continue

    # Compute ratio
    ratio = max( [fi] + neigh_flows ) / min( [fi] + neigh_flows )

    if ratio > FLOW_JUMP_FACTOR:
        # Physiologically impossible → smooth
        fixed_flow[i] = float(sum(neigh_flows) / len(neigh_flows))

# Replace with smoothed version
flow_array = fixed_flow

# Recompute power and resistance using smoothed flow
for i in valid_path_points:
    f = flow_array[i]
    if abs(f) > 1e-8:
        power_array[i]      = pressure_np[i] * f
        resistance_array[i] = pressure_np[i] / f
    else:
        power_array[i]      = 0.0
        resistance_array[i] = 0.0

# -------------------------------------------------------------------------
# Attach Arrays
# -------------------------------------------------------------------------
def add_array(name, np_array, num_comps=1):
    vtk_arr = numpy_support.numpy_to_vtk(np_array, deep=1)
    vtk_arr.SetName(name)
    vtk_arr.SetNumberOfComponents(num_comps)
    output.GetPointData().AddArray(vtk_arr)

add_array("Flow", flow_array)
add_array("Power", power_array)
add_array("Ohm's Law Resistance", resistance_array)
add_array("Slice_Normal", normals_array, 3)
"""
pf5.RequestInformationScript = ''
pf5.RequestUpdateExtentScript = ''
pf5.PythonPath = ''
pf5.UpdatePipeline()


# -------------------------
# 8) Group the result of PF5 with the original VTU
# -------------------------
grp_pf6 = GroupDatasets(Input=[pf5, XMLUnstructuredGridReader1])
grp_pf6.UpdatePipeline()



# -------------------------
# 9) Programmable Filter 6 on grouped data
# Calculates resistance based on the Poiseuille law equation
# -------------------------
pf6 = ProgrammableFilter(Input=[grp_pf6])
pf6.Script = r"""

import vtk
import math
from collections import deque

# ----------------------------
# Configuration
# ----------------------------
CLIP_THRESHOLD_FACTOR = 1.4   # if slice max distance > factor * centerline radius => clip interior
BLOOD_VISCOSITY = 0.04      # Poise (≈ 4.0 cP, the Simvascular cgs value) -- ensure unit consistency with geometry
VESSEL_LENGTH   = 0.01        # cm
MIN_RADIUS_EPS  = 1e-6        # clamp very small radii to avoid blow-ups in resistance

VERBOSE = False
def log(msg):
    if VERBOSE:
        print(msg)

# ----------------------------
# Helpers (pure VTK/Python — no NumPy reductions)
# ----------------------------
def build_connectivity(polydata):
    # Undirected adjacency (point id → neighbor ids) from polyline cells. 
    npts = polydata.GetNumberOfPoints()
    conn = [[] for _ in range(npts)]
    for cid in range(polydata.GetNumberOfCells()):
        cell = polydata.GetCell(cid)
        m = cell.GetNumberOfPoints()
        for j in range(m - 1):
            p1 = cell.GetPointId(j)
            p2 = cell.GetPointId(j + 1)
            if p2 not in conn[p1]:
                conn[p1].append(p2)
            if p1 not in conn[p2]:
                conn[p2].append(p1)
    return conn

def bfs_path(connectivity, start, end):
    # Return one shortest path (list of point ids) from start to end; [] if none.
    if start == end:
        return [start]
    visited = set([start])
    parent = {start: -1}
    dq = deque([start])
    while dq:
        u = dq.popleft()
        for nb in connectivity[u]:
            if nb not in visited:
                visited.add(nb)
                parent[nb] = u
                if nb == end:
                    # reconstruct path
                    path = [nb]
                    while parent[path[-1]] != -1:
                        path.append(parent[path[-1]])
                    path.reverse()
                    return path
                dq.append(nb)
    return []

def normalize(vx, vy, vz):
    n = math.sqrt(vx*vx + vy*vy + vz*vz)
    if n > 0.0 and math.isfinite(n):
        return (vx/n, vy/n, vz/n)
    return (1.0, 0.0, 0.0)  # safe fallback

def average_and_max_edge_radius(points_vtk, center_tuple):

    # Compute (average distance, max distance) from center to polyline vertices (edge points).
    # No NumPy reductions.
    
    if points_vtk is None:
        return (0.0, 0.0)
    n = points_vtk.GetNumberOfPoints()
    if n == 0:
        return (0.0, 0.0)
    cx, cy, cz = center_tuple
    total = 0.0
    maxd  = 0.0
    pt = [0.0, 0.0, 0.0]
    for i in range(n):
        points_vtk.GetPoint(i, pt)
        dx = pt[0] - cx
        dy = pt[1] - cy
        dz = pt[2] - cz
        d = math.sqrt(dx*dx + dy*dy + dz*dz)
        total += d
        if d > maxd:
            maxd = d
    return (total / float(n), maxd)

# ----------------------------
# Inputs
# ----------------------------
grouped_input = self.GetInputDataObject(0, 0)
if grouped_input is None:
    raise RuntimeError("Expected a MultiBlock input at port 0.")

centerline = grouped_input.GetBlock(0)
surface    = grouped_input.GetBlock(1)

if centerline is None or centerline.GetPoints() is None:
    raise RuntimeError("Centerline (block 0) missing or has no points.")
if surface is None:
    raise RuntimeError("Surface (block 1) missing.")

num_points = centerline.GetNumberOfPoints()
log(f"[Info] Centerline points: {num_points}")

pd = centerline.GetPointData()
filtered_radius_vtk = pd.GetArray("Filtered_Radius")
blanking_vtk        = pd.GetArray("blanking")
is_inlet_vtk        = pd.GetArray("IsInlet")
is_outlet_vtk       = pd.GetArray("IsOutlet")

if filtered_radius_vtk is None or blanking_vtk is None:
    raise RuntimeError("Required arrays not found: 'Filtered_Radius' and/or 'blanking'.")

# ----------------------------
# Use IsInlet / IsOutlet to build ordered paths
# ----------------------------
connectivity = build_connectivity(centerline)

inlet_ids  = []
outlet_ids = []
if is_inlet_vtk is not None and is_outlet_vtk is not None:
    for i in range(num_points):
        if int(is_inlet_vtk.GetTuple1(i)) == 1:
            inlet_ids.append(i)
        if int(is_outlet_vtk.GetTuple1(i)) == 1:
            outlet_ids.append(i)
else:
    # If either array is missing, process all points with a trivial order.
    inlet_ids  = []
    outlet_ids = []

paths = []
if inlet_ids and outlet_ids:
    for s in inlet_ids:
        for t in outlet_ids:
            path = bfs_path(connectivity, s, t)
            if path:
                paths.append(path)
    if not paths:
        # Fallback if flags present but no path found
        paths = [list(range(num_points))]
        log("[Warn] Inlet/Outlet flags present but no path found; using index order.")
else:
    # No flags: single trivial "path" using index order
    paths = [list(range(num_points))]
    log("[Warn] 'IsInlet'/'IsOutlet' missing; using index order for normals.")

# ----------------------------
# Compute per-point normals from ordered paths
# ----------------------------
path_normal = {}     # pid -> (nx, ny, nz)
path_weight = {}     # pid -> path length (to prefer normals from longer paths)
path_order  = {}     # pid -> order index along its chosen path
branch_id   = {}     # pid -> outlet id (point id) of the path this pid belongs to

pt = [0.0, 0.0, 0.0]
qt = [0.0, 0.0, 0.0]

for path in paths:
    plen = len(path)
    if plen == 0:
        continue

    # Determine outlet id for this path (last index if flags exist)
    outlet_for_path = path[-1] if outlet_ids else -1

    for k, pid in enumerate(path):
        centerline.GetPoint(pid, pt)
        if k < plen - 1:
            pid_next = path[k + 1]
            centerline.GetPoint(pid_next, qt)
            vx, vy, vz = (qt[0] - pt[0]), (qt[1] - pt[1]), (qt[2] - pt[2])
        elif plen >= 2:
            pid_prev = path[k - 1]
            centerline.GetPoint(pid_prev, qt)
            vx, vy, vz = (pt[0] - qt[0]), (pt[1] - qt[1]), (pt[2] - qt[2])
        else:
            # isolated single-point path
            vx, vy, vz = (1.0, 0.0, 0.0)

        nx, ny, nz = normalize(vx, vy, vz)

        # Only overwrite if this comes from a longer (or equal longer) path
        if (pid not in path_normal) or (plen > path_weight.get(pid, -1)):
            path_normal[pid] = (nx, ny, nz)
            path_weight[pid] = plen
            path_order[pid]  = k
            branch_id[pid]   = outlet_for_path

# Set of points to process = union of all paths
valid_path_points = set()
for pth in paths:
    for pid in pth:
        valid_path_points.add(pid)

# ----------------------------
# Prepare output arrays
# ----------------------------
avg_arr = vtk.vtkDoubleArray()
avg_arr.SetName("Slice_Average_Radius")
avg_arr.SetNumberOfComponents(1)
avg_arr.SetNumberOfTuples(num_points)
avg_arr.FillComponent(0, 0.0)

clip_arr = vtk.vtkIntArray()
clip_arr.SetName("Use_Clipping_Flag")
clip_arr.SetNumberOfComponents(1)
clip_arr.SetNumberOfTuples(num_points)
clip_arr.FillComponent(0, 0)

res_arr = vtk.vtkDoubleArray()
res_arr.SetName("Poiseuille_Resistance")
res_arr.SetNumberOfComponents(1)
res_arr.SetNumberOfTuples(num_points)
res_arr.FillComponent(0, 0.0)

# Optional: diagnostic/utility arrays
order_arr = vtk.vtkIntArray()
order_arr.SetName("PathOrder")
order_arr.SetNumberOfComponents(1)
order_arr.SetNumberOfTuples(num_points)
order_arr.FillComponent(0, -1)

branch_arr = vtk.vtkIntArray()
branch_arr.SetName("BranchOutletId")  # store outlet pid this point flows to
branch_arr.SetNumberOfComponents(1)
branch_arr.SetNumberOfTuples(num_points)
branch_arr.FillComponent(0, -1)

on_path_flag = vtk.vtkIntArray()
on_path_flag.SetName("OnAnyPathFlag")
on_path_flag.SetNumberOfComponents(1)
on_path_flag.SetNumberOfTuples(num_points)
on_path_flag.FillComponent(0, 0)

# ----------------------------
# Main loop over path points
# ----------------------------
for pid in sorted(valid_path_points):
    on_path_flag.SetTuple1(pid, 1)
    if pid in path_order:
        order_arr.SetTuple1(pid, int(path_order[pid]))
    if pid in branch_id:
        branch_arr.SetTuple1(pid, int(branch_id[pid]))

    # Center point
    p = [0.0, 0.0, 0.0]
    centerline.GetPoint(pid, p)
    p_tuple = (float(p[0]), float(p[1]), float(p[2]))

    # Path-based normal (forward along path)
    nrm = path_normal.get(pid, (1.0, 0.0, 0.0))

    # Centerline radius (for clipping heuristic)
    r_centerline = float(filtered_radius_vtk.GetTuple1(pid))
    if not (r_centerline > 0.0) or not math.isfinite(r_centerline):
        r_centerline = 1.0

    # Slice plane
    plane = vtk.vtkPlane()
    plane.SetOrigin(p_tuple)
    plane.SetNormal(nrm)

    cutter = vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputData(surface)
    cutter.Update()
    cut_output = cutter.GetOutput()

    # Keep closest contour loop to center point
    region = vtk.vtkPolyDataConnectivityFilter()
    region.SetInputData(cut_output)
    region.SetExtractionModeToClosestPointRegion()
    region.SetClosestPoint(p_tuple)
    region.Update()
    filtered_slice = region.GetOutput()

    pts = filtered_slice.GetPoints()
    avg_r = 0.0

    if pts is not None and pts.GetNumberOfPoints() > 0:
        avg_unclipped, max_unclipped = average_and_max_edge_radius(pts, p_tuple)

        use_clip = (int(blanking_vtk.GetTuple1(pid)) != 0) or (max_unclipped > CLIP_THRESHOLD_FACTOR * r_centerline)
        clip_arr.SetTuple1(pid, 1 if use_clip else 0)

        if use_clip:
            sphere = vtk.vtkSphere()
            sphere.SetCenter(p_tuple)
            sphere.SetRadius(r_centerline)

            clipper = vtk.vtkClipPolyData()
            clipper.SetInputData(filtered_slice)
            clipper.SetClipFunction(sphere)
            clipper.InsideOutOn()  # keep inside sphere
            clipper.Update()
            clipped = clipper.GetOutput()

            cpts = clipped.GetPoints()
            if cpts is not None and cpts.GetNumberOfPoints() > 0:
                avg_clipped, _ = average_and_max_edge_radius(cpts, p_tuple)
                avg_r = avg_clipped
            else:
                avg_r = avg_unclipped
        else:
            avg_r = avg_unclipped
    else:
        clip_arr.SetTuple1(pid, 0)
        avg_r = 0.0

    # Store average radius
    avg_arr.SetTuple1(pid, float(avg_r))

    # Poiseuille Resistance: R = 8 μ L / (π r^4)  (Unit-consistent!)
    r_safe = avg_r if avg_r > MIN_RADIUS_EPS else MIN_RADIUS_EPS
    R = (8.0 * BLOOD_VISCOSITY * VESSEL_LENGTH) / (math.pi * (r_safe ** 4))
    if avg_r <= MIN_RADIUS_EPS:
        R = 0.0
    res_arr.SetTuple1(pid, float(R))

# ----------------------------
# Attach arrays & output
# ----------------------------
pd.AddArray(avg_arr)
pd.AddArray(clip_arr)
pd.AddArray(res_arr)
pd.AddArray(order_arr)
pd.AddArray(branch_arr)
pd.AddArray(on_path_flag)

output.ShallowCopy(centerline)

"""
pf6.OutputDataSetType = 'vtkPolyData'
pf6.RequestInformationScript = ''
pf6.RequestUpdateExtentScript = ''
pf6.PythonPath = ''
pf6.UpdatePipeline()

# -------------------------
# 10) Programmable Filter 7
# Calculates Resistance metrics per line segment
# -------------------------
pf7 = ProgrammableFilter(Input=[pf6])
pf7.Script = r"""
import vtk
import math

# --------------------------------
# Configuration
# --------------------------------
# Radius source array (on centerline point-data)
RADIUS_ARRAY_NAME       = "Slice_Average_Radius"      # primary choice
RADIUS_FALLBACK_NAME    = "Filtered_Radius" # optional fallback if primary missing

POINT_SPACING_CM        = 0.01   # cm per point step (as requested)
BLOOD_VISCOSITY         = 0.04   # Poise (≈ 4 cP)
MIN_RADIUS_EPS          = 1e-6   # cm, clamp tiny radii to avoid blow-ups
FLOW_EPS                = 1e-12  # to avoid division by zero in ΔP/mean(Flow)

VERBOSE = False
def log(msg):
    if VERBOSE:
        print(msg)

# --------------------------------
# Connectivity and segment utilities (no NumPy)
# --------------------------------
def build_connectivity(polydata):
     # Undirected adjacency (point id -> neighbor ids) from polyline cells. 
    npts = polydata.GetNumberOfPoints()
    conn = [[] for _ in range(npts)]
    for cid in range(polydata.GetNumberOfCells()):
        cell = polydata.GetCell(cid)
        m = cell.GetNumberOfPoints()
        for j in range(m - 1):
            a = cell.GetPointId(j)
            b = cell.GetPointId(j + 1)
            if b not in conn[a]:
                conn[a].append(b)
            if a not in conn[b]:
                conn[b].append(a)
    return conn

def compute_degrees(connectivity):
    return [len(nbs) for nbs in connectivity]

def find_segments(centerline, connectivity):
  
    # Segment = maximal path between nodes with degree != 2; interior nodes have degree == 2.
    # Returns: list of segments, each a list of point IDs (including endpoints).
     
    n = centerline.GetNumberOfPoints()
    deg = compute_degrees(connectivity)
    junctions = [i for i in range(n) if deg[i] != 2]

    visited_edge = set()  # undirected edges as (min(a,b), max(a,b))
    segments = []

    def mark_edge(a, b):
        visited_edge.add((a, b) if a < b else (b, a))

    def edge_seen(a, b):
        return ((a, b) if a < b else (b, a)) in visited_edge

    # Walk outward from each junction along each neighbor, stopping at the next junction.
    for j in junctions:
        for nb in connectivity[j]:
            if edge_seen(j, nb):
                continue

            seg = [j]
            prev = j
            curr = nb
            mark_edge(prev, curr)

            while True:
                seg.append(curr)
                if deg[curr] != 2:
                    break
                nbs = connectivity[curr]
                nxt = nbs[0] if nbs[1] == prev else nbs[1]
                if edge_seen(curr, nxt):
                    break
                prev, curr = curr, nxt
                mark_edge(prev, curr)

            segments.append(seg)

    # Optional: detect/handle pure loops with all degree==2.

    return segments

def lower_quartile_mean(values):
    
    #Average of the lowest 25% of the given values.
    #Pure Python: sort ascending, take first ceil(0.25*N), average them.
    #If values is empty, return 0.0.
     
    n = len(values)
    if n == 0:
        return 0.0
    vals = list(values)
    vals.sort()
    k = int(math.ceil(0.25 * n))
    if k < 1:
        k = 1
    s = 0.0
    for i in range(k):
        s += vals[i]
    return s / float(k)

# --------------------------------
# Inputs
# --------------------------------
# Accept either a PolyData directly, or a MultiBlock with centerline in block 0
input0 = self.GetInputDataObject(0, 0)
if input0 is None:
    raise RuntimeError("No input on port 0.")

# If the input is a multiblock, grab block 0 as centerline.
if isinstance(input0, vtk.vtkMultiBlockDataSet):
    centerline = input0.GetBlock(0)
else:
    centerline = input0

if centerline is None or centerline.GetPoints() is None:
    raise RuntimeError("Centerline is missing or has no points.")

num_points = centerline.GetNumberOfPoints()
log(f"[Info] Centerline points: {num_points}")

# Arrays
pd = centerline.GetPointData()

# Radii
radius_arr = pd.GetArray(RADIUS_ARRAY_NAME)
if radius_arr is None and RADIUS_FALLBACK_NAME:
    radius_arr = pd.GetArray(RADIUS_FALLBACK_NAME)
if radius_arr is None:
    raise RuntimeError(
        f"Neither '{RADIUS_ARRAY_NAME}' nor '{RADIUS_FALLBACK_NAME}' was found on centerline point-data."
    )

# WSS (Filtered_AverageWSS)
wss_arr = pd.GetArray("Filtered_AverageWSS")
if wss_arr is None:
    log("[Warn] 'Filtered_AverageWSS' array not found. Segment_Avg_Filtered_AverageWSS will be 0.0.")

# ΔP/mean(Flow)
pressure_arr = pd.GetArray("pressure")  # case-sensitive
flow_arr     = pd.GetArray("Flow")      # case-sensitive
if pressure_arr is None or flow_arr is None:
    log("[Warn] 'pressure' and/or 'Flow' array not found. Segment_Ohm's_Resistance will be 0.")

# --------------------------------
# Build segments
# --------------------------------
connectivity = build_connectivity(centerline)
segments = find_segments(centerline, connectivity)
log(f"[Info] Found {len(segments)} segments")

# --------------------------------
# Prepare output arrays (per-point, segment-wise constant values)
# --------------------------------
seg_id_arr = vtk.vtkIntArray()
seg_id_arr.SetName("Segment_ID")
seg_id_arr.SetNumberOfComponents(1)
seg_id_arr.SetNumberOfTuples(num_points)
seg_id_arr.FillComponent(0, -1)

seg_lq25_r_arr = vtk.vtkDoubleArray()
seg_lq25_r_arr.SetName("Segment_LQ25_AvgRadius")
seg_lq25_r_arr.SetNumberOfComponents(1)
seg_lq25_r_arr.SetNumberOfTuples(num_points)
seg_lq25_r_arr.FillComponent(0, 0.0)

seg_len_arr = vtk.vtkDoubleArray()
seg_len_arr.SetName("Segment_Length_cm")
seg_len_arr.SetNumberOfComponents(1)
seg_len_arr.SetNumberOfTuples(num_points)
seg_len_arr.FillComponent(0, 0.0)

seg_res_poiseuille_arr = vtk.vtkDoubleArray()
seg_res_poiseuille_arr.SetName("Segment_Poiseuille_Resistance")
seg_res_poiseuille_arr.SetNumberOfComponents(1)
seg_res_poiseuille_arr.SetNumberOfTuples(num_points)
seg_res_poiseuille_arr.FillComponent(0, 0.0)

seg_res_datadriven_arr = vtk.vtkDoubleArray()
seg_res_datadriven_arr.SetName("Segment_Ohm's_Resistance")  # ΔP / mean(Flow)
seg_res_datadriven_arr.SetNumberOfComponents(1)
seg_res_datadriven_arr.SetNumberOfTuples(num_points)
seg_res_datadriven_arr.FillComponent(0, 0.0)

seg_wss_mean_arr = vtk.vtkDoubleArray()
seg_wss_mean_arr.SetName("Segment_Avg_Filtered_AverageWSS")
seg_wss_mean_arr.SetNumberOfComponents(1)
seg_wss_mean_arr.SetNumberOfTuples(num_points)
seg_wss_mean_arr.FillComponent(0, 0.0)

# --------------------------------
# Compute per-segment values and paint them on all points in the segment
# --------------------------------
for seg_idx, seg in enumerate(segments):
    # Collect (radius, pid) pairs for finite radii
    rpairs = []
    for pid in seg:
        r = float(radius_arr.GetTuple1(pid))
        if math.isfinite(r):
            rpairs.append((r, pid))

    # If no finite radii, everything will fall back via MIN_RADIUS_EPS later
    n_fin = len(rpairs)

    # Sort ascending by radius
    rpairs.sort(key=lambda t: t[0])

    # Determine how many are in the lowest quartile
    if n_fin > 0:
        k = int(math.ceil(0.25 * n_fin))
        if k < 1:
            k = 1
    else:
        k = 0

    # Extract the lowest-quartile subset for radius values and their PIDs
    lowest_r_vals = []
    lowest_r_pids = []
    for i in range(k):
        lowest_r_vals.append(rpairs[i][0])
        lowest_r_pids.append(rpairs[i][1])

    # 2) Lower-quartile average radius (mean of lowest_r_vals)
    if k > 0:
        lq25_avg_r = sum(lowest_r_vals) / float(k)
    else:
        lq25_avg_r = 0.0

    # 6) Segment mean of Filtered_AverageWSS across the SAME lowest-quartile radius points
    if wss_arr is not None and k > 0:
        w_sum = 0.0
        w_cnt = 0
        for pid in lowest_r_pids:
            w = float(wss_arr.GetTuple1(pid))
            if math.isfinite(w):
                w_sum += w
                w_cnt += 1
        wss_mean = (w_sum / float(w_cnt)) if w_cnt > 0 else 0.0
    else:
        wss_mean = 0.0
        # Optional: set VERBOSE=True above to see warnings
        # log(f"[Warn] seg {seg_idx}: no WSS available on LQ25 radius subset; WSS mean set to 0.")

    # 3) Segment length (as requested: number_of_points * spacing)
    seg_length_cm = float(len(seg)) * float(POINT_SPACING_CM)

    # 4) Poiseuille resistance with clamped radius
    r_safe = lq25_avg_r if lq25_avg_r > MIN_RADIUS_EPS else MIN_RADIUS_EPS
    R_poiseuille = (8.0 * BLOOD_VISCOSITY * seg_length_cm) / (math.pi * (r_safe ** 4))

    # 5) Data-driven resistance ΔP / mean(Flow)
    R_data = 0.0
    if (pressure_arr is not None) and (flow_arr is not None) and (len(seg) >= 1):
        p_start = float(pressure_arr.GetTuple1(seg[0]))
        p_end   = float(pressure_arr.GetTuple1(seg[-1]))
        q_sum = 0.0
        q_cnt = 0
        for pid in seg:
            q = float(flow_arr.GetTuple1(pid))
            if math.isfinite(q):
                q_sum += q
                q_cnt += 1
        q_mean = (q_sum / float(q_cnt)) if q_cnt > 0 else 0.0
        if abs(q_mean) > FLOW_EPS and math.isfinite(p_start) and math.isfinite(p_end):
            R_data = (p_start - p_end) / q_mean
        else:
            R_data = 0.0

    # 7) Paint values to all points in the segment
    for pid in seg:
        seg_id_arr.SetTuple1(pid, int(seg_idx))
        seg_lq25_r_arr.SetTuple1(pid, float(lq25_avg_r))
        seg_len_arr.SetTuple1(pid, float(seg_length_cm))
        seg_res_poiseuille_arr.SetTuple1(pid, float(R_poiseuille))
        seg_res_datadriven_arr.SetTuple1(pid, float(R_data))
        seg_wss_mean_arr.SetTuple1(pid, float(wss_mean))

# --------------------------------
# Attach arrays & output
# --------------------------------
pd.AddArray(seg_id_arr)
pd.AddArray(seg_lq25_r_arr)
pd.AddArray(seg_len_arr)
pd.AddArray(seg_res_poiseuille_arr)
pd.AddArray(seg_res_datadriven_arr)
pd.AddArray(seg_wss_mean_arr)

output.ShallowCopy(centerline)

"""
pf7.OutputDataSetType = 'vtkPolyData'
pf7.RequestInformationScript = ''
pf7.RequestUpdateExtentScript = ''
pf7.PythonPath = ''
pf7.UpdatePipeline()

# -------------------------
# 11) Programmable Filter 8
# Determines which 'side' all points are aligned on.
# Note, sidedness may not be geometrically accurate, as it simply groups points by which
# of the LPA or RPA they come off of (and assigns null to the MPA)
# -------------------------
pf8 = ProgrammableFilter(Input=[pf7])
pf8.Script = r"""
import vtk
from vtkmodules.util import numpy_support
import numpy as np

# Get input and output
input = self.GetInputDataObject(0, 0)
output = self.GetOutput()
output.DeepCopy(input)

num_points = input.GetNumberOfPoints()

# --- Get required array ---
generation_vtk = input.GetPointData().GetArray("Generation")

if generation_vtk is None:
    raise RuntimeError("Missing required array: 'Generation'.")

generation = numpy_support.vtk_to_numpy(generation_vtk)

# --- Build connectivity graph ---
connectivity = {i: [] for i in range(num_points)}

for cell_id in range(input.GetNumberOfCells()):
    cell = input.GetCell(cell_id)
    for j in range(cell.GetNumberOfPoints() - 1):
        p1 = cell.GetPointId(j)
        p2 = cell.GetPointId(j + 1)
        connectivity[p1].append(p2)
        connectivity[p2].append(p1)

# --- Identify trunk and branch seeds ---
trunk_points = np.where(generation < 0.1)[0]
gen1_candidates = np.where((generation > 0.9) & (generation < 1.1))[0]

branch_seeds = []
for pid in gen1_candidates:
    for neighbor in connectivity[pid]:
        if generation[neighbor] < 0.1:
            branch_seeds.append(pid)
            break

if len(branch_seeds) < 2:
    raise RuntimeError(f"Expected at least two Generation 1 branch seeds, but found {len(branch_seeds)}.")

branch1, branch2 = branch_seeds[:2]

# --- DFS restricted to one side ---
def dfs_with_tag(start, tag_label, temp_tags):
    stack = [start]
    visited = set()

    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)

        if generation[current] > 0.9 and temp_tags[current] is None:
            temp_tags[current] = tag_label

        for neighbor in connectivity[current]:
            # Do not cross trunk
            if generation[neighbor] < 0.1:
                continue
            if neighbor not in visited:
                stack.append(neighbor)

# Initialize tags
temp_tags = [None] * num_points

# Tag branches
dfs_with_tag(branch1, "Left", temp_tags)
dfs_with_tag(branch2, "Right", temp_tags)

# --- Create sidedness array ---
sidedness_vtk = vtk.vtkStringArray()
sidedness_vtk.SetName("Sidedness")
sidedness_vtk.SetNumberOfValues(num_points)

for i in range(num_points):
    if temp_tags[i] is None:
        sidedness_vtk.SetValue(i, "null")
    else:
        sidedness_vtk.SetValue(i, temp_tags[i])

output.GetPointData().AddArray(sidedness_vtk)
"""
pf8.OutputDataSetType = 'vtkPolyData'
pf8.RequestInformationScript = ''
pf8.RequestUpdateExtentScript = ''
pf8.PythonPath = ''
pf8.UpdatePipeline()

# -------------------------
# 12) Programmable Filter 9
# Renames arrays
# -------------------------
pf9 = ProgrammableFilter(Input=[pf8])
pf9.Script = r"""

from vtkmodules.util import numpy_support

# Get input
input_data = self.GetInput()

# Define array name mappings: {"old_name": "new_name"}
rename_map = {
"RegionId": "regionId",
"AverageWSS": "averageWSS",
"Blanking": "blanking",
"CapPoints": "capPoints",
"CenterlineIds": "centerlineIds",
"CrossSectionalArea": "crossSectionalArea",
"EdgeArray": "edgeArray",
"EdgePCoordArray": "edgePCoordArray",
"GlobalElementID": "globalElementID",
"GlobalNodeID": "globalNodeID",
"GroupIds": "groupIds",
"MaximumInscribedSphereRadius": "sphereRadius",
"TractIds": "tractIds",
"PathOrder": "pathOrder",
"Slice_Normal": "slice_Normal",
"Segment_ID": "segment_ID",
"Segment_LQ25_AvgRadius": "segment_avgradius",
"Segment_Data_Resistance": "z-resistance_segment_information",
"Use_Clipping_Flag": "z-clipping_Flag",
"Segment_LQ25_AvgRadius": "segment_LowerQuartile_AvgRadius",
"Segment_LQ25_AvgRadius_Initial": "z-initial_avgRadius",
"Segment_Radius_Fix_Code": "z-fix_segment_radius",
"Segment_ID": "segment_ID",
"OnAnyPathFlag": "onAnyPathFlag",
"BranchOutletId" : "z-extra_Outlet_ID",
"Power" : "power_byPressureValuePerPoint"
}

# Loop through each mapping
for old_name, new_name in rename_map.items():
    old_array = input_data.GetPointData().GetArray(old_name)
    if old_array is None:
        print(f"⚠️ Array '{old_name}' not found in PointData. Skipping.")
        continue

    # Convert to NumPy and back to VTK with new name
    array_np = numpy_support.vtk_to_numpy(old_array)
    new_array_vtk = numpy_support.numpy_to_vtk(array_np)
    new_array_vtk.SetName(new_name)

    # Add new array and remove old one
    input_data.GetPointData().AddArray(new_array_vtk)
    input_data.GetPointData().RemoveArray(old_name)
    print(f"✅ Renamed '{old_name}' to '{new_name}'.")

# Output
output.ShallowCopy(input_data)

"""
pf9.OutputDataSetType = 'vtkPolyData'
pf9.RequestInformationScript = ''
pf9.RequestUpdateExtentScript = ''
pf9.PythonPath = ''
pf9.UpdatePipeline()



# -------------------------
# 34) Programmable Filter 10
# Computes power by pressure drop
# -------------------------
pf10 = ProgrammableFilter(Input=[pf9])
pf10.Script = r"""
import vtk
from vtkmodules.util import numpy_support
import numpy as np
from collections import deque
import math

# --------------------------------
# Parameters / toggles
# --------------------------------
USE_BLANKING = False          # include all points; do not gate by Blanking
VALID_BLANKING_VALUE = 0      # unused when USE_BLANKING is False

EPS = 1e-14
FFR_MIN_BASELINE = 1e-12
UNDEFINED_FILL_VALUE = 1.0  # for undefined FFR cases (easy-to-spot sentinel)
LOG_UNDEFINED = True

# --------------------------------
# Input (PolyData or MultiBlock[0])
# --------------------------------
inp = self.GetInputDataObject(0, 0)
if inp is None:
    raise RuntimeError("No input on port 0.")

poly = inp.GetBlock(0) if isinstance(inp, vtk.vtkMultiBlockDataSet) else inp
if poly is None or poly.GetPoints() is None:
    raise RuntimeError("Input centerline PolyData missing or has no points.")

num_points = poly.GetNumberOfPoints()
pd = poly.GetPointData()

# --------------------------------
# Required arrays
# --------------------------------
is_inlet_vtk   = pd.GetArray("IsInlet")
is_outlet_vtk  = pd.GetArray("IsOutlet")
pressure_vtk   = pd.GetArray("pressure")
segment_vtk    = pd.GetArray("segment_ID")  # kept for your existing segment arrays

if is_inlet_vtk is None:
    raise RuntimeError("Missing required array: 'IsInlet'.")
if is_outlet_vtk is None:
    raise RuntimeError("Missing required array: 'IsOutlet'.")
if pressure_vtk is None:
    raise RuntimeError("Missing required array: 'pressure'.")
if segment_vtk is None:
    raise RuntimeError("Missing required array: 'segment_ID'.")

# Optional arrays (for Flow fallback)
flow_vtk          = pd.GetArray("Flow")
filtered_area_vtk = pd.GetArray("Filtered_Area")
velocity_vtk      = pd.GetArray("velocity")
blanking_vtk      = pd.GetArray("Blanking") if USE_BLANKING else None

# --------------------------------
# NumPy views
# --------------------------------
points_np   = numpy_support.vtk_to_numpy(poly.GetPoints().GetData()).astype(float, copy=False)
is_inlet    = numpy_support.vtk_to_numpy(is_inlet_vtk).astype(int, copy=False)
is_outlet   = numpy_support.vtk_to_numpy(is_outlet_vtk).astype(int, copy=False)
pressure    = numpy_support.vtk_to_numpy(pressure_vtk).astype(float, copy=False)
segment_id  = numpy_support.vtk_to_numpy(segment_vtk).astype(np.int64, copy=False)

flow = numpy_support.vtk_to_numpy(flow_vtk).astype(float, copy=False) if flow_vtk is not None else None
filtered_area = numpy_support.vtk_to_numpy(filtered_area_vtk).astype(float, copy=False) if filtered_area_vtk is not None else None
velocity      = numpy_support.vtk_to_numpy(velocity_vtk).astype(float, copy=False) if velocity_vtk is not None else None
blanking      = numpy_support.vtk_to_numpy(blanking_vtk) if blanking_vtk is not None else None

# --------------------------------
# Helpers
# --------------------------------
def add_point_array(out_poly, name, arr, num_comps=1):
    vtk_arr = numpy_support.numpy_to_vtk(arr, deep=1)
    vtk_arr.SetName(name)
    vtk_arr.SetNumberOfComponents(num_comps)
    out_poly.GetPointData().AddArray(vtk_arr)

def add_field_string_array(out_poly, name, messages):
    if not messages:
        return
    sa = vtk.vtkStringArray()
    sa.SetName(name)
    for m in messages:
        sa.InsertNextValue(m)
    out_poly.GetFieldData().AddArray(sa)

def build_connectivity(polydata):
    N = polydata.GetNumberOfPoints()
    adj = {i: [] for i in range(N)}
    for cid in range(polydata.GetNumberOfCells()):
        cell = polydata.GetCell(cid)
        npts = cell.GetNumberOfPoints()
        for j in range(npts - 1):
            a = cell.GetPointId(j)
            b = cell.GetPointId(j + 1)
            if b not in adj[a]:
                adj[a].append(b)
            if a not in adj[b]:
                adj[b].append(a)
    return adj

def multi_source_bfs(starts, adj, mask=None, with_parent=False):
    N = len(adj)
    visited = np.zeros(N, dtype=bool)
    parent = np.full(N, -1, dtype=np.int64) if with_parent else None
    dist = np.full(N, np.inf, dtype=float) if with_parent else None
    q = deque()
    for s in starts:
        if 0 <= s < N:
            if mask is None or mask[s]:
                visited[s] = True
                if with_parent:
                    dist[s] = 0.0
                    parent[s] = -1
                q.append(s)
    while q:
        u = q.popleft()
        for v in adj[u]:
            if mask is not None and not mask[v]:
                continue
            if not visited[v]:
                visited[v] = True
                if with_parent:
                    parent[v] = u
                    dist[v] = dist[u] + 1.0
                q.append(v)
    if with_parent:
        return visited, parent, dist
    else:
        return visited

# --------------------------------
# Valid subgraph & directed tree (parents/children)
# --------------------------------
adj = build_connectivity(poly)

inlet_idx  = np.where(is_inlet == 1)[0]
outlet_idx = np.where(is_outlet == 1)[0]

if inlet_idx.size == 0 or outlet_idx.size == 0:
    out = self.GetOutput()
    out.DeepCopy(poly)
    add_point_array(out, "Power_PressureDrop", np.zeros(num_points))
    add_point_array(out, "Segment_Power_PressureDrop", np.zeros(num_points))
    add_point_array(out, "Segment_PressureAverage", np.zeros(num_points))
    add_point_array(out, "Segment_FFR_perBranch", np.zeros(num_points))
    add_point_array(out, "Segment_FFR_byPressureAverage", np.zeros(num_points))
    add_point_array(out, "Parent_Branch", np.full(num_points, -1, dtype=np.int64))
    # New branch arrays
    add_point_array(out, "z-branch_ID", np.full(num_points, -1, dtype=np.int64))
    add_point_array(out, "z-branch_Parent", np.full(num_points, -1, dtype=np.int64))
    add_point_array(out, "z-branch_FFR", np.full(num_points, UNDEFINED_FILL_VALUE, dtype=float))
    raise RuntimeError("No inlets or outlets present. Arrays written as defaults.")

reach_in  = multi_source_bfs(inlet_idx,  adj)
reach_out = multi_source_bfs(outlet_idx, adj)
valid_mask = reach_in & reach_out
valid_idx = np.where(valid_mask)[0]
valid_set = set(valid_idx.tolist())

# Directed inlet-rooted tree on valid subgraph
_, parent, dist = multi_source_bfs(inlet_idx, adj, mask=valid_mask, with_parent=True)

children = {i: [] for i in range(num_points)}
for v in valid_idx:
    u = parent[v]
    if u >= 0:
        children[u].append(v)

# --------------------------------
# Flow fallback (only if no Flow is present)
# --------------------------------
tangent = np.zeros((num_points, 3), dtype=float)
need_flow_fallback = (flow is None) and (filtered_area is not None) and (velocity is not None)

if need_flow_fallback:
    for i in valid_idx:
        neigh = [n for n in adj[i] if n in valid_set]
        if not neigh:
            continue
        vecs = np.array([points_np[n] - points_np[i] for n in neigh], dtype=float)
        t = vecs.sum(axis=0)
        nrm = np.linalg.norm(t)
        if nrm > EPS and np.isfinite(nrm):
            tangent[i] = t / nrm
    flow = np.zeros(num_points, dtype=float)
    for i in valid_idx:
        if filtered_area is None or velocity is None:
            continue
        if not np.isfinite(filtered_area[i]) or filtered_area[i] <= 0.0:
            continue
        if USE_BLANKING and blanking is not None and blanking[i] != VALID_BLANKING_VALUE:
            continue
        t = tangent[i]
        if np.linalg.norm(t) <= EPS:
            continue
        v = velocity[i]
        if not np.all(np.isfinite(v)):
            continue
        flow[i] = float(np.dot(v, t)) * float(filtered_area[i])
elif flow is None:
    flow = np.zeros(num_points, dtype=float)

# --------------------------------
# Per-point: Power_PressureDrop (kept as before)
# --------------------------------
power_pd = np.zeros(num_points, dtype=float)
for i in range(num_points):
    if is_inlet[i] == 1:
        power_pd[i] = 1.0
        continue
    if not valid_mask[i]:
        continue
    up = parent[i]
    if up < 0 or not np.isfinite(pressure[i]) or not np.isfinite(pressure[up]):
        continue
    dp_i = float(pressure[up] - pressure[i])
    q_i  = float(flow[i])
    if np.isfinite(q_i):
        power_pd[i] = dp_i * q_i

# --------------------------------
# Segment-level arrays (kept for compatibility)
# --------------------------------
unique_segments = np.unique(segment_id)
seg_to_indices = {int(s): [] for s in unique_segments}
for i in range(num_points):
    seg_to_indices[int(segment_id[i])].append(i)

segment_power_pd = np.zeros(num_points, dtype=float)
segment_ffr_per_branch = np.zeros(num_points, dtype=float)
segment_pressure_avg = np.full(num_points, np.nan, dtype=float)

seg_mean_pressure = {}
seg_start_idx = {}
seg_end_idx = {}
seg_start_dist = {}

for s in unique_segments:
    inds_all = seg_to_indices[int(s)]
    if not inds_all:
        seg_mean_pressure[int(s)] = np.nan
        continue

    pvals = np.array([pressure[i] for i in inds_all], dtype=float)
    finite_mask = np.isfinite(pvals)
    seg_mean = float(np.mean(pvals[finite_mask])) if np.any(finite_mask) else np.nan
    seg_mean_pressure[int(s)] = seg_mean
    segment_pressure_avg[np.array(inds_all, dtype=int)] = seg_mean

    inds_valid = [i for i in inds_all if valid_mask[i] and np.isfinite(dist[i])]
    if not inds_valid:
        continue
    dvals = np.array([dist[i] for i in inds_valid], dtype=float)
    start_i = inds_valid[int(np.argmin(dvals))]
    end_i   = inds_valid[int(np.argmax(dvals))]

    seg_start_idx[int(s)] = start_i
    seg_end_idx[int(s)]   = end_i
    seg_start_dist[int(s)] = float(dist[start_i])

    p_start = pressure[start_i] if np.isfinite(pressure[start_i]) else np.nan
    p_end   = pressure[end_i]   if np.isfinite(pressure[end_i])   else np.nan
    q_start = flow[start_i]     if np.isfinite(flow[start_i])     else np.nan

    if np.isnan(p_start) or np.isnan(p_end) or np.isnan(q_start):
        val_ffr   = 0.0
        val_power = 0.0
    else:
        val_ffr   = (p_end / p_start) if abs(p_start) > FFR_MIN_BASELINE else 0.0
        val_power = (p_start - p_end) * q_start

    seg_inds_all = np.array(inds_all, dtype=int)
    segment_ffr_per_branch[seg_inds_all] = val_ffr
    segment_power_pd[seg_inds_all]       = val_power

# --------------------------------
# TOPOLOGICAL BRANCHING (authoritative)
# Build true branches from the inlet-rooted tree:
#   Branch = maximal chain from (root or junction) → (junction or leaf)
# --------------------------------
is_endpoint = np.zeros(num_points, dtype=bool)
for i in valid_idx:
    deg_children = len(children[i])
    if parent[i] < 0 or deg_children != 1:
        is_endpoint[i] = True

branches = []  # list of dicts with: {'id', 'up_node', 'down_node', 'nodes'}
end_node_to_branch = {}  # map downstream endpoint → branch_id

branch_id_counter = 0

# Start at all upstream endpoints; launch a branch along each outgoing child
for u in valid_idx:
    if not is_endpoint[u]:
        continue
    # For a root or junction, explore each outgoing child as a new branch
    for c in children[u] if len(children[u]) > 0 else []:
        # Walk down until next endpoint
        nodes = []
        curr = c
        while True:
            nodes.append(curr)
            if is_endpoint[curr]:
                break
            # exactly one child by construction here
            nxts = children[curr]
            if len(nxts) != 1:
                # Safety: if label got weird, break to avoid loop
                break
            curr = nxts[0]
        b = {
            'id': branch_id_counter,
            'up_node': u,
            'down_node': curr,
            'nodes': nodes
        }
        branches.append(b)
        end_node_to_branch[curr] = branch_id_counter
        branch_id_counter += 1

# Handle the special case: a root endpoint with **no children** (isolated root)
# (rare on centerlines, but just in case — make a 0-length branch)
for u in valid_idx:
    if is_endpoint[u] and len(children[u]) == 0:
        b = {
            'id': branch_id_counter,
            'up_node': u,
            'down_node': u,
            'nodes': [u]
        }
        branches.append(b)
        end_node_to_branch[u] = branch_id_counter
        branch_id_counter += 1

num_branches = len(branches)

# --------------------------------
# Build branch parent relationships
# Parent branch of branch B is the branch whose down_node equals B.up_node
# --------------------------------
branch_parent = np.full(num_branches, -1, dtype=np.int64)
upnode_to_branch = {}  # for quick lookup if needed
for b in branches:
    upnode_to_branch[b['up_node']] = b['id']

# Create mapping from endpoint nodes to branch ids (already have end_node_to_branch)
for b in branches:
    up = b['up_node']
    if up in end_node_to_branch:
        branch_parent[b['id']] = end_node_to_branch[up]  # parent branch id

# --------------------------------
# Branch mean pressures
# --------------------------------
branch_mean_pressure = np.full(num_branches, np.nan, dtype=float)
for b in branches:
    ps = np.array([pressure[i] for i in b['nodes']], dtype=float)
    m = np.isfinite(ps)
    branch_mean_pressure[b['id']] = float(np.mean(ps[m])) if np.any(m) else np.nan

# --------------------------------
# Branch FFR (z-branch_FFR): child_mean / parent_mean; trunks = 1.0; undefined = 100
# --------------------------------
z_branch_ffr_vals = np.full(num_branches, UNDEFINED_FILL_VALUE, dtype=float)
messages = []

# Identify trunk branches: those with branch_parent == -1
trunks = [b['id'] for b in branches if branch_parent[b['id']] < 0]

for bid in trunks:
    z_branch_ffr_vals[bid] = 1.0  # trunk FFR

for b in branches:
    bid = b['id']
    par = int(branch_parent[bid])
    if par < 0:
        continue  # already set to 1.0 for trunk
    child_mean = branch_mean_pressure[bid]
    parent_mean = branch_mean_pressure[par]
    if (not np.isfinite(child_mean)) or (not np.isfinite(parent_mean)) or (abs(parent_mean) <= FFR_MIN_BASELINE):
        z_branch_ffr_vals[bid] = UNDEFINED_FILL_VALUE
        if LOG_UNDEFINED:
            reason = []
            if not np.isfinite(child_mean):
                reason.append("child branch mean pressure not finite")
            if not np.isfinite(parent_mean):
                reason.append("parent branch mean pressure not finite")
            elif abs(parent_mean) <= FFR_MIN_BASELINE:
                reason.append(f"parent branch mean pressure too small (|p| ≤ {FFR_MIN_BASELINE})")
            messages.append(f"[z-branch_FFR] Branch {bid} set to {UNDEFINED_FILL_VALUE} ({'; '.join(reason)}). Parent branch {par}.")
    else:
        z_branch_ffr_vals[bid] = float(child_mean / parent_mean)

# --------------------------------
# Paint branch arrays to points (default for non-valid points)
# --------------------------------
z_branch_id_arr     = np.full(num_points, -1, dtype=np.int64)
z_branch_parent_arr = np.full(num_points, -1, dtype=np.int64)
z_branch_ffr_arr    = np.full(num_points, UNDEFINED_FILL_VALUE, dtype=float)

for b in branches:
    pts = np.array(b['nodes'], dtype=int)
    z_branch_id_arr[pts]     = b['id']
    z_branch_parent_arr[pts] = int(branch_parent[b['id']])
    z_branch_ffr_arr[pts]    = z_branch_ffr_vals[b['id']]

# --------------------------------
# OPTIONAL: keep your earlier segment-by-average arrays (for continuity)
# We will compute a simple edge-based parent for segments to fill these; note that
# the authoritative topology is given by z-branch_* arrays above.
# --------------------------------
# Edge-based (parent[v] -> v) segment transitions; choose earliest (min dist[v])
seg_parent_map   = {}
seg_boundary_dist = {}
for v in valid_idx:
    u = parent[v]
    if u < 0:
        continue
    su = int(segment_id[u]); sv = int(segment_id[v])
    if su == sv:
        continue
    dv = float(dist[v])
    if (sv not in seg_parent_map) or (dv < seg_boundary_dist[sv]):
        seg_parent_map[sv]   = su
        seg_boundary_dist[sv] = dv

segment_ffr_by_pressure_avg = np.full(num_points, UNDEFINED_FILL_VALUE, dtype=float)
parent_branch_seg_arr       = np.full(num_points, -1, dtype=np.int64)

# Determine segment roots (contain inlet points or no parent mapping)
seg_roots = set(int(segment_id[i]) for i in inlet_idx if valid_mask[int(i)])
if not seg_roots:
    seg_roots = set(int(s) for s in unique_segments if int(s) not in seg_parent_map)

# Build segment adjacency
seg_children = {}
for child_seg, par_seg in seg_parent_map.items():
    seg_children.setdefault(par_seg, []).append(child_seg)

# Traverse segment graph from roots
visited_segments = set()
queue = list(seg_roots)
for r in seg_roots:
    visited_segments.add(int(r))
    inds = np.array(seg_to_indices.get(int(r), []), dtype=int)
    if inds.size > 0:
        segment_ffr_by_pressure_avg[inds] = 1.0
        parent_branch_seg_arr[inds] = -1

while queue:
    par = int(queue.pop(0))
    for ch in seg_children.get(par, []):
        ch = int(ch)
        if ch in visited_segments:
            continue
        mean_child  = seg_mean_pressure.get(ch, np.nan)
        mean_parent = seg_mean_pressure.get(par, np.nan)
        if (not np.isfinite(mean_child)) or (not np.isfinite(mean_parent)) or (abs(mean_parent) <= FFR_MIN_BASELINE):
            val = UNDEFINED_FILL_VALUE
            if LOG_UNDEFINED:
                reason = []
                if not np.isfinite(mean_child):
                    reason.append("child segment mean pressure not finite")
                if not np.isfinite(mean_parent):
                    reason.append("parent segment mean pressure not finite")
                elif abs(mean_parent) <= FFR_MIN_BASELINE:
                    reason.append(f"parent segment mean pressure too small (|p| ≤ {FFR_MIN_BASELINE})")
                messages.append(f"[Segment_FFR_byPressureAverage] Segment {ch} set to {UNDEFINED_FILL_VALUE} ({'; '.join(reason)}). Parent segment {par}.")
        else:
            val = float(mean_child / mean_parent)
        inds = np.array(seg_to_indices.get(ch, []), dtype=int)
        if inds.size > 0:
            segment_ffr_by_pressure_avg[inds] = val
            parent_branch_seg_arr[inds] = par
        visited_segments.add(ch)
        queue.append(ch)

# --------------------------------
# Output
# --------------------------------
out = self.GetOutput()
out.DeepCopy(poly)

# Existing arrays (unchanged)
add_point_array(out, "Power_PressureDrop", power_pd)
add_point_array(out, "Segment_Power_PressureDrop", segment_power_pd)
add_point_array(out, "Segment_PressureAverage", segment_pressure_avg)
add_point_array(out, "Segment_FFR_perBranch", segment_ffr_per_branch)


# NEW: Topology-true branch arrays (authoritative)
add_point_array(out, "z-branch_ID", z_branch_id_arr)
add_point_array(out, "z-branch_Parent", z_branch_parent_arr)
add_point_array(out, "Segment_FFR_byPressureAverage", z_branch_ffr_arr)

# Log messages (undefined cases etc.)
if LOG_UNDEFINED and len(messages) > 0:
    try:
        print("\n".join(messages))
    except Exception:
        pass
    add_field_string_array(out, "BranchAndSegment_FFR_Messages", messages)
"""
pf10.OutputDataSetType = 'vtkPolyData'
pf10.RequestInformationScript = ''
pf10.RequestUpdateExtentScript = ''
pf10.PythonPath = ''
pf10.UpdatePipeline()