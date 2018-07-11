# TH2M; Small deformations, linear poroelastic, Twophase flow PP-approach

### With monolithic scheme
AddTest(
    NAME TH2M_M_cube_1e0
    PATH TH2M/M
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB cube_1e0_pcs_0_ts_*.vtu displacement displacement 1e-15 1e-15
    GLOB cube_1e0_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 1e-15
    GLOB cube_1e0_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 1e-15
    GLOB cube_1e0_pcs_0_ts_*.vtu temperature temperature 1e-15 1e-15
    GLOB cube_1e0_pcs_0_ts_*.vtu epsilon epsilon 1e-15 1e-15
    GLOB cube_1e0_pcs_0_ts_*.vtu sigma_eff sigma_eff 1e-15 1e-15
)
