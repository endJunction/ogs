# Mechanics; Small deformations Non-Local, Ehlers-damage Nooru-Mohamed (NMED)
AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NMED_Ehlers_Damage_NooruMohamed
    PATH Mechanics/EhlersDamage/Nooru_Mohamed
    EXECUTABLE ogs
    EXECUTABLE_ARGS nooru_m.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 2e-11 RELTOL 1e-13
    DIFF_DATA
    b_f__pcs_0_ts_30_t_0.030000.vtu b_f__pcs_0_ts_30_t_0.030000.vtu displacement displacement
    b_f__pcs_0_ts_30_t_0.030000.vtu b_f__pcs_0_ts_30_t_0.030000.vtu damage damage
    b_f__pcs_0_ts_30_t_0.030000.vtu b_f__pcs_0_ts_30_t_0.030000.vtu sigma_xx sigma_xx
    b_f__pcs_0_ts_30_t_0.030000.vtu b_f__pcs_0_ts_30_t_0.030000.vtu sigma_yy sigma_yy
    b_f__pcs_0_ts_30_t_0.030000.vtu b_f__pcs_0_ts_30_t_0.030000.vtu sigma_zz sigma_zz
    b_f__pcs_0_ts_30_t_0.030000.vtu b_f__pcs_0_ts_30_t_0.030000.vtu sigma_xy sigma_xy
    b_f__pcs_0_ts_40_t_0.040000.vtu b_f__pcs_0_ts_40_t_0.040000.vtu displacement displacement
    b_f__pcs_0_ts_40_t_0.040000.vtu b_f__pcs_0_ts_40_t_0.040000.vtu damage damage
    b_f__pcs_0_ts_40_t_0.040000.vtu b_f__pcs_0_ts_40_t_0.040000.vtu sigma_xx sigma_xx
    b_f__pcs_0_ts_40_t_0.040000.vtu b_f__pcs_0_ts_40_t_0.040000.vtu sigma_yy sigma_yy
    b_f__pcs_0_ts_40_t_0.040000.vtu b_f__pcs_0_ts_40_t_0.040000.vtu sigma_zz sigma_zz
    b_f__pcs_0_ts_40_t_0.040000.vtu b_f__pcs_0_ts_40_t_0.040000.vtu sigma_xy sigma_xy
    b_f__pcs_0_ts_50_t_0.050000.vtu b_f__pcs_0_ts_50_t_0.050000.vtu displacement displacement
    b_f__pcs_0_ts_50_t_0.050000.vtu b_f__pcs_0_ts_50_t_0.050000.vtu damage damage
    b_f__pcs_0_ts_50_t_0.050000.vtu b_f__pcs_0_ts_50_t_0.050000.vtu sigma_xx sigma_xx
    b_f__pcs_0_ts_50_t_0.050000.vtu b_f__pcs_0_ts_50_t_0.050000.vtu sigma_yy sigma_yy
    b_f__pcs_0_ts_50_t_0.050000.vtu b_f__pcs_0_ts_50_t_0.050000.vtu sigma_zz sigma_zz
    b_f__pcs_0_ts_50_t_0.050000.vtu b_f__pcs_0_ts_50_t_0.050000.vtu sigma_xy sigma_xy
)
