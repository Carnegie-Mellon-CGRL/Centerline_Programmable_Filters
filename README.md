The script WALT_Run performs the following operations in paraview to apply WSS, Area, Length, tortuosity, and other arrays to the centerlines:
 - Resample the centerlines with the .vtu (Applies pressure values from simulation result)
 - Filter 1 to both the resampled centerlines and the original centerlines: “Combine_Centerlines” (Gives point and cell arrays back to centerline)
 - Apply a Clean Filter at 1e-5 (Reorders points and connects some segments)
 - Filter 2: “Connect_all_lines” (Connects remaining line segments into continuous geometry)
 - Filter 3: “Centerline_Find_Segment_Lengths” (Gets length, inlet, and outlets)
 - Filter 4: “Identify_Bifurcation_Points” (Applies CapPoints, where we suspect branching regions to begin)
 - Group dataset with .vtu
 - Filter 5: “Find_Radius_at_Cap_Points” (Determines radius to apply across branching regions)
 - Cell_to_Point_Data Filter (Converts all arrays to point data arrays)
 - Group with .vtp 
 - Filter 6: “VTP_WSS_only” (Applies radius-filtered average WSS)
 - Group with .vtu again.
 - Filter 7: “VTU_Area_only” (Applies radius-filtered area)

If you want to resample the data points and apply additional arrays (Generations, Resistance, Sidedness, etc), then apply the script Resample_Run to the results of WALT_Run:
 - Filter 8: “Resample_and_Interpolate_Points” (Ensures all points are 0.01 cm apart)
 - Clean and check connectivity (Again reorders points and connects all segments)
 - Filter 9: “Generations” (Determines the generation of each vessel segment) 
 - Filter 10: “Round_up_Blanking_Values” (Rounds various values to avoid interpolation after resampling)
 - Filter 11: “inlet_outlet_and_lengths” (Determines fresh inlet, outlets, and length)
 - Calculator to get Pressure (mmHg)
 - Filter 12: “Flow_and_Resistance” (Calculates Flow, Power, and Ohm’s Law Resistance)
 - Group with vtu and apply Filter 13: “Poiseuille Resistance” (Calculates Poiseuille Law Resistance)
 - Filter 14: “Segment-wise Poiseuille Resistance (Calculates segment-wise Poiseuille Resistance)
 - Filter 15: “Sidedness” (Splits centerlines into null for the first segment, then right or left randomly)
 - Filter 16: “Rename_arrays” (Done by pre-established conventions)
 - Filter 17: “Power_by_Pressure_Drop” (Re-calculates Power based on a pressure drop, finds two versions of FFR)
