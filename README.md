The script WALT_Run performs the following operations in paraview to apply WSS, Area, Length, tortuosity, and other arrays to the centerlines:
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

If you want to resample the data points and apply additional arrays (Generations, Resistance, Sidedness, etc), then apply the script Resample_Run to the results of WALT_Run:
 - Apply Filter 7: “Resample_and_Interpolate_Points”
 - Clean and check connectivity
 - Apply Filter 8: “Generations”
 - Apply Filter 9: “Round_up_Blanking_Values”
 - Apply Filter 10: “inlet_outlet_and_lengths”
 - Apply Calculator to get Pressure (mmHg)
 - Apply Filter 11: “Flow_Power_and_Resistance”
 - Group with vtu and apply Filter 12: “True_Resistance”
 - Apply Filter 13: “Segment-wise Poiseuille Resistance
 - Apply Filter 14: “Sidedness”
 - Apply Filter 15: “Rename_arrays”
