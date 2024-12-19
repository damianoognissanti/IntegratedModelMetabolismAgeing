mutable struct CellCycleModel
    cell_cycle_sys::ODESystem
    cell_cycle_prob::ODEProblem
    p_names::Array
    u0_names::Array
    p::Array
    u0::Array
    tspan::Tuple
    porig::Array
    u0orig::Array
    tspanorig::Tuple
end

function initialiseCellCycleModel()

   porig = [
       5.029219,
       99.8,
       0.1,
       40.050851,
       2.5,
       3.0,
       15.028,
       100.0,
       1000.0,
       2.392008424,
       1.0,
       0.05,
       36.0,
       2.0105,
       50.0,
       0.0005,
       1.596,
       0.1,
       0.5,
       1.0,
       60.0,
       22.025,
       1.0,
       1.0093466795,
       0.5,
       1.0e-5,
       3.0,
       0.1,
       0.5,
       0.05,
       1.0,
       0.22,
       0.05,
       20.029,
       0.005,
       1.0,
       0.1,
       0.0001,
       100.1,
       6.0,
       4.058984,
       0.01405,
       0.368352,
       10.0,
       2.99896964,
       0.025,
       0.00700465,
       0.0033,
       0.2,
       0.43,
       600.0,
       0.013,
       0.2527,
       30.0,
       0.1,
       94288.25,
       15.0,
       0.01,
       3.169628,
       0.187825,
       1.001,
       1.001,
       1.87453125,
       0.0995,
       1.0,
       3.0,
       0.1,
       5000.0,
       1.5115,
       15.7842816,
       6.0,
       0.02,
       10.0,
       3.0,
       0.5,
       0.4032023,
       0.1,
       40.0,
       4.0,
       0.5028,
       20.404,
       0.01,
       0.30555,
       30.0,
       350.267663,
       0.5005,
       3.0,
       121.568756607,
       0.14066,
       1.0,
       2.0,
       1.0,
       28.0496125,
       85.0316195,
       0.604,
       6.0,
       21.0,
       0.16533,
       0.01,
       10.0,
       671.037047,
       0.025,
       0.214,
       0.56044,
       0.1,
       0.5,
       1.0,
       0.2267,
       24.99,
       0.013,
       0.01,
       50.0
    ]

   u0orig = [
       0.0,
       70.2499994605,
       0.0,
       0.0055095267392,
       0.000515873455915,
       1.10259641365,
       0.0,
       1.13255806709,
       0.0,
       285.604158558,
       0.0616005253745,
       0.327003199617,
       277.885559446,
       0.0268170651527,
       0.89919603461,
       338.4307188605,
       0.0410087523757,
       349.937883622,
       52.41210074,
       0.390588672104,
       620.374730215,
       3421.56928094,
       253.1708654155,
       0.0,
       198.280009914,
       22.84051778295,
       0.0,
       1.622362084625,
       571.275477765,
       0.0,
       47.7610290701,
       16.5140084901,
       263.70599151,
       107.0,
       0.612700299515
    ]

   p_names = [
       "ki_Cdc14",
       "kpp_Swi5_Clb2",
       "kd_Sic1p",
       "v0_Mcm1",
       "Kp_Mcm1",
       "ka_APC_Cdc14",
       "Kp_Sic1",
       "kpp_APC_Clb5",
       "Kp_Swi5",
       "kp_Swi5",
       "kcd_Clb3_Sic1",
       "kdp_Far1p",
       "kpp_Cln2_Far1p",
       "kpp_Cln2_Whi5",
       "kdp_Swi5_Cdc14",
       "kd_Cln3_Far1p",
       "kpp_Swe1_Hls1",
       "kd_Cln2_Far1p",
       "kcd_Cln3_Far1p",
       "Kp_Clb3",
       "kp_Mcm1",
       "Kp_Cln2",
       "Kp_Clb5",
       "kpp_APC_Clb2",
       "kd_APC_p",
       "kpp_Cdc14_MEN",
       "kcf_SBF_Whi5",
       "kd_Swe1_p",
       "kcd_Cln2_Far1p",
       "kd_Clb3_APC",
       "kI_Clb5_Hog1",
       "kd_APC",
       "kdd_Far1p",
       "kI_Cln2_Hog1",
       "kd_Clb3",
       "kcd_Clb5_Sic1",
       "kd_Sic1",
       "kI_Clb2_Hog1",
       "kI_Swe1_Hog1",
       "n1",
       "kpp_SBF_Clb2",
       "kp_Cln3",
       "kp_Clb3",
       "kpp_Far1",
       "kp_Far1",
       "kd_Clb2_p",
       "kpp_Cdc14_Clb2",
       "kdd_Far1",
       "kdp_Clb5_Sic1_Hp",
       "kpp_Cln2_Sic1",
       "kcf_Clb5_Sic1_Hp",
       "kd_Swi5",
       "kd_MBF",
       "kcf_Clb3_Sic1",
       "kd_Clb5_Sic1",
       "Kp_APC",
       "kcf_Clb2_Sic1",
       "kd_Whi5",
       "kp_Sic1",
       "kd_Clb5",
       "kd_Clb5_APC",
       "Kpp_Cln2_Whi5",
       "kp_Clb5",
       "kpp_Swe1_Clb2",
       "nutrition_factor",
       "n_Mcm1",
       "kd_Mih1",
       "Kp_Clb2",
       "kdp_Sic1_Hp",
       "kdp_Clb2",
       "kcf_Cln3_Far1p",
       "kd_Far1",
       "kpp_Clb5_Sic1_Hog1",
       "kp_MBF",
       "kdp_SBF",
       "kpp_Cln3_Whi5",
       "kd_Swe1",
       "Kpp_Cln2_Sic1",
       "n_SBF",
       "ka_Cdc14_APC",
       "V0_Mcm1",
       "kd_Clb3_Sic1",
       "kd_Mcm1",
       "kcf_Clb5_Sic1",
       "kpp_Clb2",
       "kpp_Swi5_Clb5",
       "n_Clb3",
       "kp_Clb2",
       "kp_basal_Far1",
       "kcd_Clb2_Sic1",
       "kcd_Clb5_Sic1_Hp",
       "K_MBF",
       "kp_Cln2",
       "Kpp_Cln3_Whi5",
       "kd_Clb2_APC",
       "kcf_Cln2_Far1p",
       "kpp_Clb5_Sic1",
       "kp_Whi5",
       "kd_Whi5p",
       "kpp_Sic1_Hog1",
       "kp_APC",
       "kd_Clb2",
       "kp_Mih1",
       "kp_Swe1",
       "kd_Far1p",
       "kd_Clb2_Sic1",
       "kdp_Whi5",
       "kd_Cln2",
       "Kpp_Clb5_Sic1",
       "kd_Swi5_p",
       "kd_Cln3",
       "compartment_1"
    ]

   u0_names = [
       "Sic1_Hp",
       "Cln3",
       "Fus3",
       "Clb5",
       "SBF_p",
       "Swi5_p",
       "Hog1_PP",
       "Clb2_p",
       "Far1_p",
       "Sic1",
       "SBF",
       "Clb2",
       "Clb3_Sic1",
       "MBF",
       "Clb5_Sic1",
       "Cdc14_p",
       "APC_p",
       "SBF_Whi5",
       "Whi5_p",
       "Cln2",
       "Whi5",
       "Cdc14",
       "Far1",
       "Cln3_Far1_p",
       "Mcm1",
       "Clb2_Sic1",
       "Clb5_Sic1_Hp",
       "Clb3",
       "Swi5",
       "Cln2_Far1_p",
       "APC",
       "Swe1",
       "Swe1_p",
       "Mih1",
       "Sic1_p"
    ]

    cell_cycle_sys = createCellCycleModel()
    tspanorig = (1,350)
    cell_cycle_prob = createCellCycleProblem(cell_cycle_sys, u0orig, tspanorig, porig)

    return CellCycleModel(cell_cycle_sys, cell_cycle_prob, p_names, u0_names, deepcopy(porig), deepcopy(u0orig), deepcopy(tspanorig), porig, u0orig, tspanorig)

end

function getCCParamIndex(cellcycle::CellCycleModel, name)
    return findfirst(x->x==name, cellcycle.p_names)
end

function getCCStateIndex(cellcycle::CellCycleModel, name)
    return findfirst(x->x==name, cellcycle.u0_names)
end

function getCCParamValue(cellcycle::CellCycleModel, name)
    return cellcycle.p[getCCParamIndex(cellcycle, name)]
end

function getCCStateValue(cellcycle::CellCycleModel, name)
    return cellcycle.u0[getCCStateIndex(cellcycle, name)]
end

function getCCParamOrigValue(cellcycle::CellCycleModel, name)
    return cellcycle.porig[getCCParamIndex(cellcycle, name)]
end

function getCCStateOrigValue(cellcycle::CellCycleModel, name)
    return cellcycle.u0orig[getCCStateIndex(cellcycle, name)]
end

function setCCParamValue(cellcycle::CellCycleModel, name, value)
    cellcycle.p[getCCParamIndex(cellcycle, name)] = value
    return 
end

function setCCStateValue(cellcycle::CellCycleModel, name, value)
    cellcycle.u0[getCCStateIndex(cellcycle, name)] = value
    return
end

function updateCCProblem(cellcycle::CellCycleModel)
    # Updates cell cycle ode problem with current values set
    cellcycle.cell_cycle_prob = createCellCycleProblem(cellcycle.cell_cycle_sys, cellcycle.u0, cellcycle.tspan, cellcycle.p)  
    return
end

function resetCCProblem(cellcycle::CellCycleModel)
    # Updates cell cycle ode problem with original values set
    cellcycle.u0=deepcopy(cellcycle.u0orig)
    cellcycle.p=deepcopy(cellcycle.porig)
    cellcycle.tspan=deepcopy(cellcycle.tspanorig)
    cellcycle.cell_cycle_prob = createCellCycleProblem(cellcycle.cell_cycle_sys, cellcycle.u0, cellcycle.tspan, cellcycle.p)  
    return
end
