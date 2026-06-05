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
 - Apply a Clean at lowest threshold (Again reorders points and connects all segments)
 - Filter 9: “Generations” (Determines the generation of each vessel segment)
 - Filter 10: “Round_up_Blanking_Values” (Rounds various values to avoid interpolation after resampling)
 - Filter 11: “inlet_outlet_and_lengths” (Determines fresh inlet, outlets, and length)
 - Calculator to get Pressure (mmHg)
 - Group with .vtu
 - Filter 12: “VolumetricFlowRate” (Calculates flow by the surface integral of the velocity vector)
 - Filter 13: “Flow_Power_and_Resistance” (Calculates Power, and Ohm’s Law Resistance)
 - Group with vtu and apply Filter 14: “True_Resistance” (Calculates Poiseuille Law Resistance)
 - Filter 15: “Segment Resistance (Calculates segment-wise Poiseuille Resistance and other segment-wise functions)
 - Filter 16: “Sidedness” (Splits centerlines into null for the first segment, then right or left randomly)
 - Filter 17: “Rename_arrays” (Done by pre-established conventions)
 - Filter 18: “Power_by_pressure_drop” (segment-wise power)
 - Filter 19: “Segment_Power_and_Flow” (calculates segment flow and power)
 - Filter 20: “Delete_arrays” (Removes extra arrays to declutter)
 - Filter 21: “Threshold_FFR” (Thresholds twice, once at 0.8…)
 - Filter 22: “Threshold_FFR” (...and once at 0.9)
 - Filter 23: “Threshold_Power” (Thresholds at a set value of power)
 - Filter 24: “PowerandFFROverlap” (Thresholds values within either the FFR_0.8 or FFR_0.9 arrays that coincide with the Power cutoff)


The centerline arrays are divided roughly into 'geometrical markers', 'sliced results', 'length and distance measurements', 'Hemodynamic measurements', and 'segment-wise measurements'.  Note, what you see below is NOT the order they appear in Paraview (they are alphabetized there):

Geometrical Markers
IsInlet -> Has a value of 1 for the point determined as the inlet centerline point and a value of 0 for all other points.
IsOutlet -> Has a value of 1 for points determined to be outlet centerline points and a value of 0 for all other points.
IsMedianPoint -> Has a value of 1 for points designated as the ‘midpoint’ of each segment, where the surrounding points are used for flow averaging
Generation -> Dictates the generation of the vessel.  The first branch after the inlet is 0, and every time there is a new segment it increments the generation by 1.
Sidedness -> Records all points in the RPA as "Right", all points in the LPA as "Left", and any other points (usually points in the MPA) as "null"

Sliced Results
Filtered_Area -> Slices through the full area of the 3D geometry at all points that are NOT in clipping flag criteria (regions likely to be branching areas), and at points in the clipping flag criteria it restricts the slicing area to a given radius
Filtered_Average_WSS -> Averages the sliced WSS for a given point on the surrounding vtp.  Same clipping criteria as above
Filtered_Radius -> The radius that it will restrict the slice to if within a clipping-flag point
Slice_Average_Radius -> Average Radius at given point to surrounding geometry

Length and Distance Measurements
LinearDistanceFromInlet -> Takes a distance from the point marked by IsInlet to the current centerline point.  Calculation is the distance of the linear line connecting the two points.
CenterlineLength -> Equivalent of Arc Length.  Sums the individual distances from the current point to the inlet by calculating individual distances from that point to all points that lead to the inlet.
Tortuosity -> CenterlineLength divided by LinearDistanceFromInlet

Hemodynamic Measurements
Pressure (mmHg) -> pressure sampled from simulation results onto centerlines and converted to mmHg
Flow -> Determines surface integral of each slice and takes dot product of it by the velocity
Power -> Multiplies pressure at each point by the flow at that point.
Power_PressureDrop -> The pressure drop from the prior point to the current multiplied by flow
Poiseuille_Resistance -> Vascular resistance calculated at each point based on Poiseuille formula.  Radius is the average radius at that point, length is distance from prior point to current point, and viscosity is 0.04 cP
Ohm's Law Resistance -> Vascular resistance calculated at each point based on (P2-P1)/Flow = Resistance
Threshold_FFR_0.8 -> Marks segments with a value of 1 that have a FFR less than 0.8.
Threshold_FFR_0.9 -> Marks segments with a value of 1 that have a FFR less than 0.9. 
Threshold_FFR0.9_Power_Overlap -> Marks segments that are both within the Threshold FFR < 0.9 and a manually assigned value of power threshold.
Threshold_Power_highestMulti-modelQuartile -> Marks segments with a value of 1 that are above the calculated threshold of power marking the highest quartile of all models
Threshold_Power_highestMulti-modelQuartile -> Marks segments with a value of 1 that are above the calculated threshold of power marking the highest quartile of each model
VolumetricFlowRate -> Same calculation as flow but without fixing the directionality (positive or negative based on dot product)

Segment-wise Measurements
Note, a given segment is defined as the points in between two branching regions.
Segment_AverageWSS - >Determines the average WSS of all points on a given segment and applies that average across the entire centerline
Segment_Length_cm -> Measures the distance from the start of the segment to the end of the segment
Segment_Ohm's_Resistance -> Takes the pressure at the first point of the segment minus the pressure at the last point of the segment times the average flow throughout the centerline.
Segment_Poiseuille_Resistance -> Finds the lower quartile of radius values for a given segment and the length of that segment, and calculates a segment-wise Poiseuille Resistance.
Segment_Power_PressureDrop -> Determines the pressure at the start of the segment minus the pressure at the end of the segment times the average flow of the points surrounding the point marked by IsMedian
Segment_PressureAverage -> Takes the average pressure across the entire segment in cgs units
Segment_PressureAverage (mmHg) -> Average pressure across segment in mmHg
Segment_FFR_perBranch -> Determines the ratio of the pressure at the end of the segment to the pressure at the start of the segment
Segment_FFR_byPressureAverage -> Determines the ratio of the average pressure across the current branch to the average pressure across its parent branch (edited) 
Segment_Flow_offset_outlet -> Averages the flow at a point 15 points upstream and the surrounding 5 points to either side of the last point in the segment
Segment_Power_offset_outlet -> Multiplies the pressure drop of the segment by the flow 15 points upstream of the last point in the segment.
Segment_Power_offset_outlet_average -> Multiplies the pressure drop of the segment by the average of the flow 15 points upstream of the last point in the segment and the surrounding 5 points to either side.

There are also less significant arrays that have value when coding or sorting information.  They will be described shortly here:


average_pressure -> carry over from simulation results
average_speed -> carry over from simulation results
blanking -> Simvascular-assigned branching region
capPoints -> Beginning and end points of branching regions
centerlineIds -> Paraview-assigned sorting of points.  Meaningless.
filtered_Radius_Initial -> First attempt to find radius.
globalElementID -> randomly assigns a unique ID tag to each centerline point
globalNodeID -> generated by Simvascular.  Unknown purpose
groupIDs -> groups together centerline points.  Unknown sorting
onAnyPathFlag -> tracks connectivity of centerline points to inlet
pathOrder -> organizes points based on relation to inlet (0).
power_byPressureValuesPerPoint -> calculates power as pressure at that point times flow at that point.  
pressure -> pressure in cgs units
segment_ID -> divides centerlines into segments.  Centerlines with same segment_ID belong to same segment
segment_LowerQuartile_AvgRadius -> Computes average radius at all points on a segment, finds lowest 25% of those values, and gets the average of them.
segment_Power_Average -> Old method of calculating power
slice_Normal -> the direction vector of the slice
sphereRadius -> carry over from Simvascular extraction
timeDeriv -> carry over from Simvascular extraction
tractIds -> similar to groupIds
vWSS -> leftover result of simulation resampling
velocity -> velocity at all points based on simulation results
vinplane_traction -> leftover of Simvascular extraction
vtkValidPointMask -> identifies valid points on centerline
z-branch_ID -> extra method of identifying segment ID
z-branch_parent -> Lists parent branch of each segment
z-clipping_Flag -> lists which points should be clipped with a value of 1
z-extra_FFR -> Additional calculation of FFR
z-extra_Outlet_ID -> Lists the centerline ID of a path that arrives on an outlet from current position.
