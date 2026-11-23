STANDALONE CENTERLINE POST-PROCESSING
==================================

Files
-----
- WALT_Run_standalone.py: VTP centerline + VTU simulation + (optional) VTP surface WSS.
  Produces centerline topology/length/generation, radius, clot-aware fluid
  cross-sectional pressure/flow/area, and wall-surface WSS when supplied.
- Resample_Run_standalone.py: VTP centerline + VTU simulation.
  Produces resampled centerlines, geometry fields, clot-aware pressure/flow/area,
  power/resistance, Poiseuille estimates, segment metrics, FFR, and thresholds.

Dependencies
------------
Use a normal Python environment (not ParaView/PVPython):

    python -m pip install numpy vtk

Either method of running works:

A) Spyder / normal Python IDE:
   Set CENTERLINES_PATH, VTU_PATH, and for WALT, SURFACE_VTP_PATH at the
   beginning of the appropriate script. Then click Run.

B) Windows PowerShell / Linux shell:

    python -u WALT_Run_standalone.py --centerlines "C:\\path\\centerlines.vtp" --vtu "C:\\path\\results.vtu" --surface-vtp "C:\\path\\results.vtp"

    python -u Resample_Run_standalone.py --centerlines "C:\\path\\centerlines.vtp" --vtu "C:\\path\\results.vtu"

Results
-------
Automatically saved in a directory named "centerline_results" beside
centerlines.vtp. Override using --output-dir "C:\\path\\my_results" or
OUTPUT_DIR near the top of each script. This does not modify the input files
and does not add anything to a running ParaView session.

Each completed stage writes an indexed .vtp file, e.g.:
- WALT_09_09_WSS_wall_integral.vtp
- WALT_10_10_pressure_flow_area_integrals.vtp
- Resample_07_07_pressure_flow_area_integrals.vtp
- Resample_13_13_FFR.vtp

Final files:
- centerline_results/WALT_centerlines_final.vtp
- centerline_results/Resample_centerlines_final.vtp

Per-run progress/error logs (written after every stage):
- WALT_progress.json
- Resample_progress.json

If an error occurs, completed VTP checkpoints remain available for inspection
in ParaView. A failed stage does not overwrite preceding checkpoint files.

Accuracy / performance controls
-------------------------------
The fast versions' adaptive branch-aware sampling is retained:
--section-step 0.10  (in geometry units, typically cm)
--wss-step 0.12      (in geometry units, typically cm)
--section-refinement 2
--wss-refinement 1
--full-resolution    (direct sampling at every centerline point, slower)
--no-stage-files     (write only major stage files, conserving disk space)

Caveats
-------
* Verify that geometry units are cm and VTU pressure is in dyn/cm^2; the
  embedded source code converts pressure to mmHg with /1333.22.
* Pressure is area-averaged over actual fluid cross-section triangles;
  flow is integrated normal velocity. Clot voids are excluded by construction.
* Area-effective radius and Poiseuille resistance are approximations for
  irregular/noncircular, clot-obstructed lumens.
* Sparse interpolation can miss very localized stenosis peaks; use a smaller
  --section-step or --full-resolution to compare.
* WALT wall-shear stress is only calculated when a results surface VTP
  containing a WSS array is provided.
* If the input centerline lacks explicit trustworthy Left/Right sidedness,
  Resample assigns the first two initial daughter branches in deterministic
  point-ID order. Verify anatomical Left/Right labels before analyses.
* Neither pipeline has yet been validated on the user's full patient files.
  They were exercised with synthetic VTP/VTU test data containing a void and,
  for Resample, a Y-shaped bifurcation.
