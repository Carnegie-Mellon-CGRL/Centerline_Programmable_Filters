This is the full pathway and order in which these filters are applied:
 - Resample the centerlines to get pressure
 - Apply filter 1 to both the resampled centerlines and the original centerlines:
 - Filter 1: “Combine_Centerlines”
 - Clean and check Connectivity so all are connected
 - Filter 2: “Centerline_Find_Segment_Lengths”
 - Filter 3: “Identify_Bifurcation_Points”
 - Group dataset with .vtu
 - Filter 4: “Find_Radius_at_Cap_Points”
 - Cell_to_Point_Data Filter
 - Group with .vtp 
 - Filter 5: “VTP_WSS_only”
 - Group with .vtu again.
 - Filter 6: “VTU_Area_only”
		Note: Output data type should always be set to VTKPolyData 

