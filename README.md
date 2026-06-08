# WALT: Workflow for Analysis of Lumen Transport

## Overview

**WALT** is a ParaView-based pipeline for processing vascular centerlines and integrating computational simulation outputs (e.g., pressure, velocity, wall shear stress) to derive geometric and hemodynamic metrics.

The workflow is designed to produce reproducible, point-wise and segment-wise measurements along vascular networks, including:

- Geometric descriptors (length, tortuosity, branching structure)
- Surface-derived quantities (WSS, area, radius)
- Hemodynamic metrics (pressure, flow, resistance, power)
- Functional indicators (fractional flow reserve, segment-level behavior)

The pipeline consists of two stages:

- **`WALT_Run`** — Core processing and metric extraction  
- **`Resample_Run`** — Optional refinement and higher-order analysis  

---

# Pipeline Description

## 1. `WALT_Run`: Core Processing Workflow

This pipeline maps simulation data onto centerlines and computes foundational geometric and hemodynamic quantities.

---

### 1.1 Data Mapping and Initialization

1. **Resample centerlines with `.vtu`**
   - Applies simulation-derived pressure to centerline points using **closest-point linear interpolation**

2. **Filter 1: `Combine_Centerlines`**
   - Combines resampled centerlines with original centerlines  
   - Restores all original **point and cell arrays**

3. **Clean Filter (tolerance = `1e-5`)**
   - Reorders points
   - Reconnects partially disconnected segments

4. **Filter 2: `Connect_all_lines`**
   - Connects all remaining segments into a **single continuous centerline network**

---

### 1.2 Topology and Connectivity Analysis

5. **Filter 3: `Centerline_Find_Segment_Lengths`**
   - Identifies:
     - **Inlet**: The point with minimal distance to the inflow boundary and associated with the first `CenterlineID`
     - **Outlets**: All points connected to only **one neighboring point**
   - Computes path lengths using **Breadth-First Search (BFS)**:
     - Distances are accumulated point-to-point from inlet to each outlet
     - Connectivity determines traversal order

6. **Filter 4: `Identify_Bifurcation_Points`**
   - Identifies branching regions via:
     - Points with **≥ 3 connections**
     - Segments where the `Blanking` array begins
   - Assigns **CapPoints**, representing the **start and end of bifurcation regions**

---

### 1.3 Radius and Surface Quantities

7. **Group dataset with `.vtu`**

8. **Filter 5: `Find_Radius_at_Cap_Points`**
   - Determines local vessel radius by:
     - Slicing each segment at CapPoints
     - Measuring the **maximum radial distance from centerline point to surface geometry**

9. **CellToPointData**
   - Converts all arrays to **point-associated data**

10. **Group with `.vtp`**

---

### 1.4 Surface Sampling

11. **Filter 6: `VTP_WSS_only`**
   - Computes **radius-filtered average WSS** by:
     - Slicing the surface (`.vtp`) at each centerline point
     - Averaging all **non-zero WSS values**
   - In branching regions:
     - Slice is restricted to points within the radius determined by **Filter 5**

12. **Group with `.vtu`**

13. **Filter 7: `VTU_Area_only`**
   - Computes **cross-sectional area** by:
     - Slicing volumetric data (`.vtu`)
     - Integrating area across slice
   - In branching regions:
     - Slice is restricted using previously computed radius

---

## 2. `Resample_Run`: Advanced Processing Workflow (Optional)

This step improves spatial consistency and computes higher-order physiological metrics.

---

### 2.1 Resampling and Connectivity

8. **Filter 8: `Resample_and_Interpolate_Points`**
   - Ensures uniform spacing (**0.01 cm**) along centerlines by:
     - Retaining existing points
     - Interpolating new points between them from inlet to each outlet

9. **Clean Filter (minimal threshold)**
   - Reorders and reconnects all segments post-resampling

---

### 2.2 Structural and Topological Metrics

10. **Filter 9: `Generations`**
   - Assigns generation values:
     - Starts at **0 after the inlet**
     - Increments by **1 at each branching point** (≥ 3 connections)

11. **Filter 10: `Round_up_Blanking_Values`**
   - Adjusts values to prevent interpolation artifacts introduced during resampling

12. **Filter 11: `inlet_outlet_and_lengths`**
   - Recomputes:
     - Inlet and outlet locations
     - Lengths via the same BFS method as Filter 3
   - Ensures consistency after interpolation

---

### 2.3 Hemodynamic Calculations

13. **Pressure conversion**
   - Converts mapped pressure to **mmHg**

14. **Group with `.vtu`**

15. **Filter 12: `VolumetricFlowRate`**
   - Computes flow by:
     - Performing a **surface integral of velocity**
     - Evaluated via slicing at each centerline point

16. **Filter 13: `Resistance`**
   - Computes **Ohm’s Law resistance**:
     - \( R = \frac{\Delta P}{Q} \)

17. **Filter 14: `True_Resistance`**
   - Computes **Poiseuille resistance**:
     - \( R = \frac{8 \mu L}{\pi r^4} \)
     - \( \mu = 4.0 \) cP  
     - \( L = 0.01 \) cm  
     - \( r \) = average slice radius at each point  

---

### 2.4 Segment-Based Computation

18. **Filter 15: `Segment_Resistance_WSS_and_segmentID`**
   - Defines segments:
     - A segment is the region **between two bifurcation points**
     - `segmentID` increments at each branch
   - Computes:
     - **Segment-wise Ohm’s resistance**:
       - Pressure drop from segment start to end × flow near segment end
     - **Segment-wise Poiseuille resistance**:
       - Uses total segment length and average radius
     - **Segment-average WSS**

---

### 2.5 Classification and Cleanup

19. **Filter 16: `Sidedness`**
   - Assigns:
     - `Left` and `Right` branches at first bifurcation (arbitrarily assigned)
     - `Null` for proximal (pre-branch) segments

20. **Filter 17: `Rename_arrays`**
   - Standardizes naming conventions

---

### 2.6 Power and FFR Analysis

21. **Filter 18: `Power_by_pressure_drop`**
   - Computes:
     - **Point-wise power**:
       - Pressure drop from previous point × flow at current point
     - **FFR metrics**:
       - Inter-segment and intra-segment definitions

22. **Filter 19: `Segment_Power_and_Flow`**
   - Computes:
     - Flow = average flow near **segment outlet**
     - Power = segment pressure drop × outlet flow

---

### 2.7 Thresholding and Finalization

23. **Filter 20: `Delete_arrays`**
   - Removes intermediate arrays

24. **Filters 21–24: Thresholding**
   - FFR thresholds (<0.8, <0.9)
   - Power thresholds
   - Combined FFR-power classification

---

# Data Structure and Array Definitions

## Geometrical Markers

- **`IsInlet`** — 1 at inlet point, 0 elsewhere  
- **`IsOutlet`** — 1 at outlet points (only one connection)  
- **`IsMedianPoint`** — midpoint of segment used for local averaging  
- **`Generation`** — increments at each bifurcation  
- **`Sidedness`** — Left / Right / Null classification  

---

## Sliced Results (Surface-Based)

- **`Filtered_Area`**  
  - Area computed via slicing  
  - In branching regions, restricted to radius-defined neighborhood  

- **`Filtered_Average_WSS`**  
  - Mean WSS from slice  
  - Only non-zero values included  

- **`Filtered_Radius`**  
  - Radius used to restrict slicing  

- **`Slice_Average_Radius`**  
  - Mean radius from all surface points in slice  

---

## Length and Distance Measurements

- **`LinearDistanceFromInlet`**  
  - Straight-line distance from inlet  

- **`CenterlineLength`**  
  - Path (arc) length computed via point-to-point accumulation  

- **`Tortuosity`**  
  - Ratio:
    ```
    Tortuosity = CenterlineLength / LinearDistanceFromInlet
    ```

---

## Hemodynamic Measurements

- **`Pressure (mmHg)`** — mapped simulation pressure  
- **`Flow`** — directional flow via velocity dot product  
- **`VolumetricFlowRate`** — unsigned flow (direction ignored)  

- **`Power`** — \( P \times Q \)  
- **`Power_PressureDrop`** — \( \Delta P \times Q \)  

- **`Poiseuille_Resistance`**  
  - Uses local radius and step-wise length  

- **`Ohms_Law_Resistance`**  
  - \( \frac{\Delta P}{Q} \)

---

### Threshold Indicators

- **`Threshold_FFR_0.8`** — FFR < 0.8  
- **`Threshold_FFR_0.9`** — FFR < 0.9  
- **`Threshold_FFR0.9_Power_Overlap`** — combined criteria  
- **`Threshold_Power_highestMulti-modelQuartile`**  
  - Identifies segments in highest quartile of power  

---

## Segment-wise Measurements

A **segment is defined as points between two bifurcations**.

- **`Segment_AverageWSS`**  
  - Mean WSS across all points in segment  

- **`Segment_Length_cm`**  
  - Total segment length  

- **`Segment_Ohm's_Resistance`**  
  - Pressure drop across segment × average flow  

- **`Segment_Poiseuille_Resistance`**  
  - Uses:
    - Segment length
    - Average radius (often lower quartile to reduce bias)

- **`Segment_Power_PressureDrop`**  
  - Segment pressure drop × average flow near **IsMedianPoint**

- **`Segment_FFR_perBranch`**  
  - Ratio:
    ```
    P_end / P_start
    ```

- **`Segment_FFR_byPressureAverage`**  
  - Ratio of average branch pressure to parent branch pressure  

- **`Segment_Flow_offset_outlet`**  
  - Flow averaged:
    - 15 points upstream of outlet  
    - ±5 neighboring points  

- **`Segment_Power_offset_outlet`**  
  - Segment ΔP × flow 15 points upstream  

- **`Segment_Power_offset_outlet_average`**  
  - Same as above but using averaged flow window  

---

## Auxiliary Arrays

Used for intermediate computation or bookkeeping:

- Simulation-derived:
  - `pressure`, `velocity`, `average_pressure`, `average_speed`

- Branching:
  - `blanking`, `capPoints`

- Structural:
  - `segment_ID`, `centerlineIds`, `groupIDs`, `tractIds`

- Connectivity:
  - `onAnyPathFlag`, `pathOrder`

- Misc:
  - `slice_Normal`, `vtkValidPointMask`, `sphereRadius`

- Experimental:
  - `z-*` arrays (alternate or legacy calculations)

---

## Notes

- Arrays in ParaView appear **alphabetically**, not in pipeline order  
- Segment definitions strictly follow **bifurcation topology**  
- Radius-restricted slicing is critical to avoid **branch contamination**  
- Resampling ensures **uniform spatial resolution**, improving numerical stability  

---

## Citation

If you use this workflow in academic research, please cite this repository and associated publications as appropriate.

---

## Contact

For questions or collaboration inquiries, please contact the maintainers.
