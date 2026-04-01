from paraview.simple import *
import os
import numpy as np
import pandas as pd
from vtk.util.numpy_support import vtk_to_numpy

# --------------------------------------------------------
# User input
# --------------------------------------------------------

centerlines_path = r"C:/Users/Admin/Documents/Simvascular_Files/0-Centerlines/0-Final_Results_Mar30_0.2Threshold/JL22572_Centerlines.vtp"

# Output Excel path (same name as centerline)
output_dir = os.path.dirname(centerlines_path)
base_name = os.path.splitext(os.path.basename(centerlines_path))[0]
output_csv = os.path.join(output_dir, base_name + ".csv")

# --------------------------------------------------------
# Load Centerline
# --------------------------------------------------------

centerlines_src = XMLPolyDataReader(FileName=[centerlines_path])
centerlines_src.UpdatePipeline()

# Fetch data to client (important!)
cl_data = servermanager.Fetch(centerlines_src)

point_data = cl_data.GetPointData()
num_points = cl_data.GetNumberOfPoints()

# --------------------------------------------------------
# Convert point data to pandas dataframe
# --------------------------------------------------------

df = pd.DataFrame()

df = pd.DataFrame()

for i in range(point_data.GetNumberOfArrays()):
    arr = point_data.GetArray(i)

    if arr is None:
        continue  

    name = arr.GetName()
    np_arr = vtk_to_numpy(arr)

    # Scalar array
    if np_arr.ndim == 1:
        df[name] = np_arr

    # Vector / tensor array
    else:
        for j in range(np_arr.shape[1]):
            df[f"{name}_{j}"] = np_arr[:, j]


print("Converted VTK point data to DataFrame")
print(df.head())

# --------------------------------------------------------
# Extracting Segment-wise data
# --------------------------------------------------------

Data = pd.DataFrame()

uniqueIDs = df["segment_ID"].unique()
uniqueIDs = np.sort(uniqueIDs)

for x in uniqueIDs:

    if pd.isna(x):
        print("Skipping NaN segments")
        continue

    segment_df = df.loc[df["segment_ID"] == x]
    columnID = f"segmentID_{int(x)}"
    Partial = segment_df.iloc[1]
    Data[columnID] = Partial

# Optional transpose if you want "one row per segment"
Data_out = Data.T
Data_out.index.name = "segment_ID"


# --------------------------------------------------------
# Export to excel
# --------------------------------------------------------

Data_out.to_csv(output_csv, index=True)

