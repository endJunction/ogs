AddTest(
    NAME HeatTransportBHE_3D_beier_sandbox
    PATH Parabolic/T/3D_Beier_sandbox
    EXECUTABLE ogs
    EXECUTABLE_ARGS beier_sandbox.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    beier_sandbox_pcs_0_ts_10_t_600.000000.vtu
    beier_sandbox_pcs_0_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    beier_sandbox_pcs_0_ts_10_t_600.000000.vtu beier_sandbox_pcs_0_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)
