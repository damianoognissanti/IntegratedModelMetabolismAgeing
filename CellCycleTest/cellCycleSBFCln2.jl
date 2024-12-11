ModelingToolkit.@variables t Sic1_Hp(t) Cln3(t) Fus3(t) Clb5(t) SBF_p(t) Swi5_p(t) Hog1_PP(t) Clb2_p(t) Far1_p(t) Sic1(t) SBF(t) Clb2(t) Clb3_Sic1(t) MBF(t) Clb5_Sic1(t) Cdc14_p(t) APC_p(t) SBF_Whi5(t) Whi5_p(t) Cln2(t) Whi5(t) Cdc14(t) Far1(t) Cln3_Far1_p(t) Mcm1(t) Clb2_Sic1(t) Clb5_Sic1_Hp(t) Clb3(t) Swi5(t) Cln2_Far1_p(t) APC(t) Swe1(t) Swe1_p(t) Mih1(t) Sic1_p(t)

state_array = [Sic1_Hp, Cln3, Fus3, Clb5, SBF_p, Swi5_p, Hog1_PP, Clb2_p, Far1_p, Sic1, SBF, Clb2, Clb3_Sic1, MBF, Clb5_Sic1, Cdc14_p, APC_p, SBF_Whi5, Whi5_p, Cln2, Whi5, Cdc14, Far1, Cln3_Far1_p, Mcm1, Clb2_Sic1, Clb5_Sic1_Hp, Clb3, Swi5, Cln2_Far1_p, APC, Swe1, Swe1_p, Mih1, Sic1_p]

ModelingToolkit.@parameters ki_Cdc14 kpp_Swi5_Clb2 kd_Sic1p v0_Mcm1 Kp_Mcm1 ka_APC_Cdc14 Kp_Sic1 kpp_APC_Clb5 Kp_Swi5 kp_Swi5 kcd_Clb3_Sic1 kdp_Far1p kpp_Cln2_Far1p kpp_Cln2_Whi5 kdp_Swi5_Cdc14 kd_Cln3_Far1p kpp_Swe1_Hls1 kd_Cln2_Far1p kcd_Cln3_Far1p Kp_Clb3 kp_Mcm1 Kp_Cln2 Kp_Clb5 kpp_APC_Clb2 kd_APC_p kpp_Cdc14_MEN kcf_SBF_Whi5 kd_Swe1_p kcd_Cln2_Far1p kd_Clb3_APC kI_Clb5_Hog1 kd_APC kdd_Far1p kI_Cln2_Hog1 kd_Clb3 kcd_Clb5_Sic1 kd_Sic1 kI_Clb2_Hog1 kI_Swe1_Hog1 n1 kpp_SBF_Clb2 kp_Cln3 kp_Clb3 kpp_Far1 kp_Far1 kd_Clb2_p kpp_Cdc14_Clb2 kdd_Far1 kdp_Clb5_Sic1_Hp kpp_Cln2_Sic1 kcf_Clb5_Sic1_Hp kd_Swi5 kd_MBF kcf_Clb3_Sic1 kd_Clb5_Sic1 Kp_APC kcf_Clb2_Sic1 kd_Whi5 kp_Sic1 kd_Clb5 kd_Clb5_APC Kpp_Cln2_Whi5 kp_Clb5 kpp_Swe1_Clb2 nutrition_factor n_Mcm1 kd_Mih1 Kp_Clb2 kdp_Sic1_Hp kdp_Clb2 kcf_Cln3_Far1p kd_Far1 kpp_Clb5_Sic1_Hog1 kp_MBF kdp_SBF kpp_Cln3_Whi5 kd_Swe1 Kpp_Cln2_Sic1 n_SBF ka_Cdc14_APC V0_Mcm1 kd_Clb3_Sic1 kd_Mcm1 kcf_Clb5_Sic1 kpp_Clb2 kpp_Swi5_Clb5 n_Clb3 kp_Clb2 kp_basal_Far1 kcd_Clb2_Sic1 kcd_Clb5_Sic1_Hp K_MBF kp_Cln2 Kpp_Cln3_Whi5 kd_Clb2_APC kcf_Cln2_Far1p kpp_Clb5_Sic1 kp_Whi5 kd_Whi5p kpp_Sic1_Hog1 kp_APC kd_Clb2 kp_Mih1 kp_Swe1 kd_Far1p kd_Clb2_Sic1 kdp_Whi5 kd_Cln2 Kpp_Clb5_Sic1 kd_Swi5_p kd_Cln3 compartment_1

parameter_array = [ki_Cdc14, kpp_Swi5_Clb2, kd_Sic1p, v0_Mcm1, Kp_Mcm1, ka_APC_Cdc14, Kp_Sic1, kpp_APC_Clb5, Kp_Swi5, kp_Swi5, kcd_Clb3_Sic1, kdp_Far1p, kpp_Cln2_Far1p, kpp_Cln2_Whi5, kdp_Swi5_Cdc14, kd_Cln3_Far1p, kpp_Swe1_Hls1, kd_Cln2_Far1p, kcd_Cln3_Far1p, Kp_Clb3, kp_Mcm1, Kp_Cln2, Kp_Clb5, kpp_APC_Clb2, kd_APC_p, kpp_Cdc14_MEN, kcf_SBF_Whi5, kd_Swe1_p, kcd_Cln2_Far1p, kd_Clb3_APC, kI_Clb5_Hog1, kd_APC, kdd_Far1p, kI_Cln2_Hog1, kd_Clb3, kcd_Clb5_Sic1, kd_Sic1, kI_Clb2_Hog1, kI_Swe1_Hog1, n1, kpp_SBF_Clb2, kp_Cln3, kp_Clb3, kpp_Far1, kp_Far1, kd_Clb2_p, kpp_Cdc14_Clb2, kdd_Far1, kdp_Clb5_Sic1_Hp, kpp_Cln2_Sic1, kcf_Clb5_Sic1_Hp, kd_Swi5, kd_MBF, kcf_Clb3_Sic1, kd_Clb5_Sic1, Kp_APC, kcf_Clb2_Sic1, kd_Whi5, kp_Sic1, kd_Clb5, kd_Clb5_APC, Kpp_Cln2_Whi5, kp_Clb5, kpp_Swe1_Clb2, nutrition_factor, n_Mcm1, kd_Mih1, Kp_Clb2, kdp_Sic1_Hp, kdp_Clb2, kcf_Cln3_Far1p, kd_Far1, kpp_Clb5_Sic1_Hog1, kp_MBF, kdp_SBF, kpp_Cln3_Whi5, kd_Swe1, Kpp_Cln2_Sic1, n_SBF, ka_Cdc14_APC, V0_Mcm1, kd_Clb3_Sic1, kd_Mcm1, kcf_Clb5_Sic1, kpp_Clb2, kpp_Swi5_Clb5, n_Clb3, kp_Clb2, kp_basal_Far1, kcd_Clb2_Sic1, kcd_Clb5_Sic1_Hp, K_MBF, kp_Cln2, Kpp_Cln3_Whi5, kd_Clb2_APC, kcf_Cln2_Far1p, kpp_Clb5_Sic1, kp_Whi5, kd_Whi5p, kpp_Sic1_Hog1, kp_APC, kd_Clb2, kp_Mih1, kp_Swe1, kd_Far1p, kd_Clb2_Sic1, kdp_Whi5, kd_Cln2, Kpp_Clb5_Sic1, kd_Swi5_p, kd_Cln3, compartment_1]

Eq = [
Differential(t)(Sic1_Hp) ~  1/(compartment_1) * (compartment_1*kcd_Clb5_Sic1_Hp*Clb5_Sic1_Hp - compartment_1*kdp_Sic1_Hp*Sic1_Hp - compartment_1*kcf_Clb5_Sic1_Hp*Clb5*Sic1_Hp + compartment_1*kpp_Sic1_Hog1*Hog1_PP*Sic1),
Differential(t)(Cln3) ~  1/(compartment_1) * (compartment_1*kp_Cln3 + compartment_1*kcd_Cln3_Far1p*Cln3_Far1_p - compartment_1*kd_Cln3*Cln3 - compartment_1*kcf_Cln3_Far1p*Cln3*Far1_p + compartment_1*kpp_Cln2_Far1p*Cln2*Cln3_Far1_p),
Differential(t)(Fus3) ~  1/(compartment_1) * (0),
Differential(t)(Clb5) ~  1/(compartment_1) * ((compartment_1*(kpp_Clb5_Sic1*(Kpp_Cln2_Sic1^n1)*(Clb5^n1) + kpp_Clb5_Sic1*(Clb5^n1)*(Cln2^n1) + kpp_Cln2_Sic1*(Kpp_Clb5_Sic1^n1)*(Cln2^n1) + kpp_Cln2_Sic1*(Clb5^n1)*(Cln2^n1))*Clb5_Sic1) / ((Kpp_Clb5_Sic1^n1 + Clb5^n1)*(Kpp_Cln2_Sic1^n1 + Cln2^n1)) + (compartment_1*kp_Clb5*MBF) / ((Kp_Clb5 + MBF)*(1 + kI_Clb5_Hog1*Hog1_PP)) + compartment_1*kcd_Clb5_Sic1*Clb5_Sic1 + compartment_1*kcd_Clb5_Sic1_Hp*Clb5_Sic1_Hp - compartment_1*kd_Clb5*Clb5 - compartment_1*kcf_Clb5_Sic1*Clb5*Sic1 - compartment_1*kcf_Clb5_Sic1_Hp*Clb5*Sic1_Hp - compartment_1*kd_Clb5_APC*Clb5*APC),
Differential(t)(SBF_p) ~  1/(compartment_1) * (-compartment_1*kdp_SBF*Cdc14_p*SBF_p + compartment_1*kpp_SBF_Clb2*Clb2*SBF),
Differential(t)(Swi5_p) ~  1/(compartment_1) * (-compartment_1*kd_Swi5_p*Swi5_p - compartment_1*kdp_Swi5_Cdc14*Swi5_p*Cdc14_p + compartment_1*kpp_Swi5_Clb2*Clb2*Swi5 + compartment_1*kpp_Swi5_Clb5*Clb5*Swi5),
Differential(t)(Hog1_PP) ~  1/(compartment_1) * (0),
Differential(t)(Clb2_p) ~  1/(compartment_1) * (-compartment_1*kd_Clb2_p*Clb2_p - compartment_1*kdp_Clb2*Mih1*Clb2_p + compartment_1*kpp_Clb2*Swe1*Clb2),
Differential(t)(Far1_p) ~  1/(compartment_1) * (compartment_1*kcd_Cln2_Far1p*Cln2_Far1_p + compartment_1*kcd_Cln3_Far1p*Cln3_Far1_p - compartment_1*kd_Far1p*Far1_p - compartment_1*kdp_Far1p*Far1_p - compartment_1*kcf_Cln2_Far1p*Cln2*Far1_p - compartment_1*kcf_Cln3_Far1p*Cln3*Far1_p - compartment_1*kdd_Far1p*Cln2*Far1_p + compartment_1*kpp_Far1*Far1*Fus3),
Differential(t)(Sic1) ~  1/(compartment_1) * ((-compartment_1*(kpp_Clb5_Sic1*(Kpp_Cln2_Sic1^n1)*(Clb5^n1) + kpp_Clb5_Sic1*(Clb5^n1)*(Cln2^n1) + kpp_Cln2_Sic1*(Kpp_Clb5_Sic1^n1)*(Cln2^n1) + kpp_Cln2_Sic1*(Clb5^n1)*(Cln2^n1))*Sic1) / ((Kpp_Clb5_Sic1^n1 + Clb5^n1)*(Kpp_Cln2_Sic1^n1 + Cln2^n1)) + (compartment_1*kp_Sic1*Swi5) / (Kp_Sic1 + Swi5) + compartment_1*kcd_Clb2_Sic1*Clb2_Sic1 + compartment_1*kcd_Clb3_Sic1*Clb3_Sic1 + compartment_1*kcd_Clb5_Sic1*Clb5_Sic1 - compartment_1*kd_Sic1*Sic1 + compartment_1*kdp_Sic1_Hp*Sic1_Hp - compartment_1*kcf_Clb2_Sic1*Clb2*Sic1 - compartment_1*kcf_Clb3_Sic1*Clb3*Sic1 - compartment_1*kcf_Clb5_Sic1*Clb5*Sic1 - compartment_1*kpp_Sic1_Hog1*Hog1_PP*Sic1),
Differential(t)(SBF) ~  1/(compartment_1) * (((Kpp_Cln3_Whi5*kpp_Cln2_Whi5*(Cln2^n_SBF) + kpp_Cln2_Whi5*Cln3*(Cln2^n_SBF) + kpp_Cln3_Whi5*(Kpp_Cln2_Whi5^n_SBF)*Cln3 + kpp_Cln3_Whi5*Cln3*(Cln2^n_SBF))*compartment_1*SBF_Whi5) / ((Kpp_Cln3_Whi5 + Cln3)*(Kpp_Cln2_Whi5^n_SBF + Cln2^n_SBF)) - compartment_1*kcf_SBF_Whi5*SBF*Whi5 + compartment_1*kdp_SBF*Cdc14_p*SBF_p - compartment_1*kpp_SBF_Clb2*Clb2*SBF),
Differential(t)(Clb2) ~  1/(compartment_1) * ((compartment_1*kp_Clb2*Mcm1) / ((Kp_Clb2 + Mcm1)*(1 + kI_Clb2_Hog1*Hog1_PP)) + compartment_1*kcd_Clb2_Sic1*Clb2_Sic1 - compartment_1*kd_Clb2*Clb2 - compartment_1*kcf_Clb2_Sic1*Clb2*Sic1 - compartment_1*kd_Clb2_APC*APC*Clb2 + compartment_1*kdp_Clb2*Mih1*Clb2_p - compartment_1*kpp_Clb2*Swe1*Clb2),
Differential(t)(Clb3_Sic1) ~  1/(compartment_1) * (-compartment_1*kcd_Clb3_Sic1*Clb3_Sic1 - compartment_1*kd_Clb3_Sic1*Clb3_Sic1 + compartment_1*kcf_Clb3_Sic1*Clb3*Sic1),
Differential(t)(MBF) ~  1/(compartment_1) * ((compartment_1*kp_MBF*(Cln2^n1)) / (K_MBF^n1 + Cln2^n1) - compartment_1*kd_MBF*MBF),
Differential(t)(Clb5_Sic1) ~  1/(compartment_1) * ((-compartment_1*(kpp_Clb5_Sic1*(Kpp_Cln2_Sic1^n1)*(Clb5^n1) + kpp_Clb5_Sic1*(Clb5^n1)*(Cln2^n1) + kpp_Cln2_Sic1*(Kpp_Clb5_Sic1^n1)*(Cln2^n1) + kpp_Cln2_Sic1*(Clb5^n1)*(Cln2^n1))*Clb5_Sic1) / ((Kpp_Clb5_Sic1^n1 + Clb5^n1)*(Kpp_Cln2_Sic1^n1 + Cln2^n1)) - compartment_1*kcd_Clb5_Sic1*Clb5_Sic1 - compartment_1*kd_Clb5_Sic1*Clb5_Sic1 + compartment_1*kdp_Clb5_Sic1_Hp*Clb5_Sic1_Hp + compartment_1*kcf_Clb5_Sic1*Clb5*Sic1 - compartment_1*kpp_Clb5_Sic1_Hog1*Clb5_Sic1*Hog1_PP),
Differential(t)(Cdc14_p) ~  1/(compartment_1) * (-compartment_1*ki_Cdc14*Cdc14_p + compartment_1*ka_Cdc14_APC*APC*Cdc14 + compartment_1*kpp_Cdc14_Clb2*Clb2*Cdc14 + compartment_1*kpp_Cdc14_MEN*Cdc14_p*Cdc14),
Differential(t)(APC_p) ~  1/(compartment_1) * (-compartment_1*kd_APC_p*APC_p - compartment_1*ka_APC_Cdc14*APC_p*Cdc14_p + compartment_1*kpp_APC_Clb2*APC*Clb2 + compartment_1*kpp_APC_Clb5*Clb5*APC),
Differential(t)(SBF_Whi5) ~  1/(compartment_1) * ((-(Kpp_Cln3_Whi5*kpp_Cln2_Whi5*(Cln2^n_SBF) + kpp_Cln2_Whi5*Cln3*(Cln2^n_SBF) + kpp_Cln3_Whi5*(Kpp_Cln2_Whi5^n_SBF)*Cln3 + kpp_Cln3_Whi5*Cln3*(Cln2^n_SBF))*compartment_1*SBF_Whi5) / ((Kpp_Cln3_Whi5 + Cln3)*(Kpp_Cln2_Whi5^n_SBF + Cln2^n_SBF)) + compartment_1*kcf_SBF_Whi5*SBF*Whi5),
Differential(t)(Whi5_p) ~  1/(compartment_1) * (((Kpp_Cln3_Whi5*kpp_Cln2_Whi5*(Cln2^n_SBF) + kpp_Cln2_Whi5*Cln3*(Cln2^n_SBF) + kpp_Cln3_Whi5*(Kpp_Cln2_Whi5^n_SBF)*Cln3 + kpp_Cln3_Whi5*Cln3*(Cln2^n_SBF))*compartment_1*SBF_Whi5) / ((Kpp_Cln3_Whi5 + Cln3)*(Kpp_Cln2_Whi5^n_SBF + Cln2^n_SBF)) - compartment_1*kd_Whi5p*Whi5_p - compartment_1*kdp_Whi5*Cdc14_p*Whi5_p + compartment_1*(kpp_Cln2_Whi5*Cln2 + kpp_Cln3_Whi5*Cln3)*Whi5),
Differential(t)(Cln2) ~  1/(compartment_1) * ((compartment_1*kp_Cln2*SBF) / ((Kp_Cln2 + SBF)*(1 + kI_Cln2_Hog1*Hog1_PP)) + compartment_1*kcd_Cln2_Far1p*Cln2_Far1_p - compartment_1*kd_Cln2*Cln2 - compartment_1*kcf_Cln2_Far1p*Cln2*Far1_p + compartment_1*kpp_Cln2_Far1p*Cln2*Cln2_Far1_p),
Differential(t)(Whi5) ~  1/(compartment_1) * (compartment_1*kp_Whi5 - compartment_1*kd_Whi5*Whi5 - compartment_1*kcf_SBF_Whi5*SBF*Whi5 + compartment_1*kdp_Whi5*Cdc14_p*Whi5_p - compartment_1*(kpp_Cln2_Whi5*Cln2 + kpp_Cln3_Whi5*Cln3)*Whi5),
Differential(t)(Cdc14) ~  1/(compartment_1) * (compartment_1*ki_Cdc14*Cdc14_p - compartment_1*ka_Cdc14_APC*APC*Cdc14 - compartment_1*kpp_Cdc14_Clb2*Clb2*Cdc14 - compartment_1*kpp_Cdc14_MEN*Cdc14_p*Cdc14),
Differential(t)(Far1) ~  1/(compartment_1) * (-compartment_1*kd_Far1*Far1 + compartment_1*kdp_Far1p*Far1_p + compartment_1*(kp_basal_Far1 + kp_Far1*Fus3) - compartment_1*kdd_Far1*Far1*Cln2 - compartment_1*kpp_Far1*Far1*Fus3),
Differential(t)(Cln3_Far1_p) ~  1/(compartment_1) * (-compartment_1*kcd_Cln3_Far1p*Cln3_Far1_p - compartment_1*kd_Cln3_Far1p*Cln3_Far1_p + compartment_1*kcf_Cln3_Far1p*Cln3*Far1_p - compartment_1*kpp_Cln2_Far1p*Cln2*Cln3_Far1_p),
Differential(t)(Mcm1) ~  1/(compartment_1) * ((compartment_1*kp_Mcm1*(Clb2^n_Mcm1)) / (Kp_Mcm1^n_Mcm1 + Clb2^n_Mcm1) + (compartment_1*v0_Mcm1*(Clb3^n_Mcm1)) / (V0_Mcm1^n_Mcm1 + Clb3^n_Mcm1) - compartment_1*kd_Mcm1*Mcm1),
Differential(t)(Clb2_Sic1) ~  1/(compartment_1) * (-compartment_1*kcd_Clb2_Sic1*Clb2_Sic1 - compartment_1*kd_Clb2_Sic1*Clb2_Sic1 + compartment_1*kcf_Clb2_Sic1*Clb2*Sic1),
Differential(t)(Clb5_Sic1_Hp) ~  1/(compartment_1) * (-compartment_1*kcd_Clb5_Sic1_Hp*Clb5_Sic1_Hp - compartment_1*kdp_Clb5_Sic1_Hp*Clb5_Sic1_Hp + compartment_1*kcf_Clb5_Sic1_Hp*Clb5*Sic1_Hp + compartment_1*kpp_Clb5_Sic1_Hog1*Clb5_Sic1*Hog1_PP),
Differential(t)(Clb3) ~  1/(compartment_1) * ((compartment_1*kp_Clb3*(Clb5^n_Clb3)) / (Kp_Clb3^n_Clb3 + Clb5^n_Clb3) + compartment_1*kcd_Clb3_Sic1*Clb3_Sic1 - compartment_1*kd_Clb3*Clb3 - compartment_1*kcf_Clb3_Sic1*Clb3*Sic1 - compartment_1*kd_Clb3_APC*APC*Clb3),
Differential(t)(Swi5) ~  1/(compartment_1) * ((compartment_1*kp_Swi5*Mcm1) / (Kp_Swi5 + Mcm1) - compartment_1*kd_Swi5*Swi5 + compartment_1*kdp_Swi5_Cdc14*Swi5_p*Cdc14_p - compartment_1*kpp_Swi5_Clb2*Clb2*Swi5 - compartment_1*kpp_Swi5_Clb5*Clb5*Swi5),
Differential(t)(Cln2_Far1_p) ~  1/(compartment_1) * (-compartment_1*kcd_Cln2_Far1p*Cln2_Far1_p - compartment_1*kd_Cln2_Far1p*Cln2_Far1_p + compartment_1*kcf_Cln2_Far1p*Cln2*Far1_p - compartment_1*kpp_Cln2_Far1p*Cln2*Cln2_Far1_p),
Differential(t)(APC) ~  1/(compartment_1) * ((compartment_1*kp_APC*Mcm1) / (Kp_APC + Mcm1) - compartment_1*kd_APC*APC + compartment_1*ka_APC_Cdc14*APC_p*Cdc14_p - compartment_1*kpp_APC_Clb2*APC*Clb2 - compartment_1*kpp_APC_Clb5*Clb5*APC),
Differential(t)(Swe1) ~  1/(compartment_1) * ((-compartment_1*kpp_Swe1_Hls1*Swe1) / (1 + kI_Swe1_Hog1*Hog1_PP) + compartment_1*kp_Swe1 - compartment_1*kd_Swe1*Swe1 - compartment_1*kpp_Swe1_Clb2*Swe1*Clb2),
Differential(t)(Swe1_p) ~  1/(compartment_1) * ((compartment_1*kpp_Swe1_Hls1*Swe1) / (1 + kI_Swe1_Hog1*Hog1_PP) - compartment_1*kd_Swe1_p*Swe1_p + compartment_1*kpp_Swe1_Clb2*Swe1*Clb2),
Differential(t)(Mih1) ~  1/(compartment_1) * (compartment_1*kp_Mih1 - compartment_1*kd_Mih1*Mih1),
Differential(t)(Sic1_p) ~  1/(compartment_1) * ((compartment_1*(kpp_Clb5_Sic1*(Kpp_Cln2_Sic1^n1)*(Clb5^n1) + kpp_Clb5_Sic1*(Clb5^n1)*(Cln2^n1) + kpp_Cln2_Sic1*(Kpp_Clb5_Sic1^n1)*(Cln2^n1) + kpp_Cln2_Sic1*(Clb5^n1)*(Cln2^n1))*Sic1) / ((Kpp_Clb5_Sic1^n1 + Clb5^n1)*(Kpp_Cln2_Sic1^n1 + Cln2^n1)) + (compartment_1*(kpp_Clb5_Sic1*(Kpp_Cln2_Sic1^n1)*(Clb5^n1) + kpp_Clb5_Sic1*(Clb5^n1)*(Cln2^n1) + kpp_Cln2_Sic1*(Kpp_Clb5_Sic1^n1)*(Cln2^n1) + kpp_Cln2_Sic1*(Clb5^n1)*(Cln2^n1))*Clb5_Sic1) / ((Kpp_Clb5_Sic1^n1 + Clb5^n1)*(Kpp_Cln2_Sic1^n1 + Cln2^n1)) - compartment_1*kd_Sic1p*Sic1_p)
]

@named osys = ODESystem(Eq, t, state_array, parameter_array)

p = [
ki_Cdc14 => 5.029219,
kpp_Swi5_Clb2 => 99.8,
kd_Sic1p => 0.1,
v0_Mcm1 => 40.050851,
Kp_Mcm1 => 2.5,
ka_APC_Cdc14 => 3.0,
Kp_Sic1 => 15.028,
kpp_APC_Clb5 => 100.0,
Kp_Swi5 => 1000.0,
kp_Swi5 => 2.392008424,
kcd_Clb3_Sic1 => 1.0,
kdp_Far1p => 0.05,
kpp_Cln2_Far1p => 36.0,
kpp_Cln2_Whi5 => 2.0105,
kdp_Swi5_Cdc14 => 50.0,
kd_Cln3_Far1p => 0.0005,
kpp_Swe1_Hls1 => 1.596,
kd_Cln2_Far1p => 0.1,
kcd_Cln3_Far1p => 0.5,
Kp_Clb3 => 1.0,
kp_Mcm1 => 60.0,
Kp_Cln2 => 22.025,
Kp_Clb5 => 1.0,
kpp_APC_Clb2 => 1.0093466795,
kd_APC_p => 0.5,
kpp_Cdc14_MEN => 1.0e-5,
kcf_SBF_Whi5 => 3.0,
kd_Swe1_p => 0.1,
kcd_Cln2_Far1p => 0.5,
kd_Clb3_APC => 0.05,
kI_Clb5_Hog1 => 1.0,
kd_APC => 0.22,
kdd_Far1p => 0.05,
kI_Cln2_Hog1 => 20.029,
kd_Clb3 => 0.005,
kcd_Clb5_Sic1 => 1.0,
kd_Sic1 => 0.1,
kI_Clb2_Hog1 => 0.0001,
kI_Swe1_Hog1 => 100.1,
n1 => 6.0,
kpp_SBF_Clb2 => 4.058984,
kp_Cln3 => 0.01405,
kp_Clb3 => 0.368352,
kpp_Far1 => 10.0,
kp_Far1 => 2.99896964,
kd_Clb2_p => 0.025,
kpp_Cdc14_Clb2 => 0.00700465,
kdd_Far1 => 0.0033,
kdp_Clb5_Sic1_Hp => 0.2,
kpp_Cln2_Sic1 => 0.43,
kcf_Clb5_Sic1_Hp => 600.0,
kd_Swi5 => 0.013,
kd_MBF => 0.2527,
kcf_Clb3_Sic1 => 30.0,
kd_Clb5_Sic1 => 0.1,
Kp_APC => 94288.25,
kcf_Clb2_Sic1 => 15.0,
kd_Whi5 => 0.01,
kp_Sic1 => 3.169628,
kd_Clb5 => 0.187825,
kd_Clb5_APC => 1.001,
Kpp_Cln2_Whi5 => 1.001,
kp_Clb5 => 1.87453125,
kpp_Swe1_Clb2 => 0.0995,
nutrition_factor => 1.0,
n_Mcm1 => 3.0,
kd_Mih1 => 0.1,
Kp_Clb2 => 5000.0,
kdp_Sic1_Hp => 1.5115,
kdp_Clb2 => 15.7842816,
kcf_Cln3_Far1p => 6.0,
kd_Far1 => 0.02,
kpp_Clb5_Sic1_Hog1 => 10.0,
kp_MBF => 3.0,
kdp_SBF => 0.5,
kpp_Cln3_Whi5 => 0.4032023,
kd_Swe1 => 0.1,
Kpp_Cln2_Sic1 => 40.0,
n_SBF => 4.0,
ka_Cdc14_APC => 0.5028,
V0_Mcm1 => 20.404,
kd_Clb3_Sic1 => 0.01,
kd_Mcm1 => 0.30555,
kcf_Clb5_Sic1 => 30.0,
kpp_Clb2 => 350.267663,
kpp_Swi5_Clb5 => 0.5005,
n_Clb3 => 3.0,
kp_Clb2 => 121.568756607,
kp_basal_Far1 => 0.14066,
kcd_Clb2_Sic1 => 1.0,
kcd_Clb5_Sic1_Hp => 2.0,
K_MBF => 1.0,
kp_Cln2 => 28.0496125,
Kpp_Cln3_Whi5 => 85.0316195,
kd_Clb2_APC => 0.604,
kcf_Cln2_Far1p => 6.0,
kpp_Clb5_Sic1 => 21.0,
kp_Whi5 => 0.16533,
kd_Whi5p => 0.01,
kpp_Sic1_Hog1 => 10.0,
kp_APC => 671.037047,
kd_Clb2 => 0.025,
kp_Mih1 => 0.214,
kp_Swe1 => 0.56044,
kd_Far1p => 0.1,
kd_Clb2_Sic1 => 0.5,
kdp_Whi5 => 1.0,
kd_Cln2 => 0.2267,
Kpp_Clb5_Sic1 => 24.99,
kd_Swi5_p => 0.013,
kd_Cln3 => 0.01,
compartment_1 => 50.0
]

u0 = [
Sic1_Hp =>  1/(compartment_1) * (0.0),
Cln3 =>  1/(compartment_1) * (70.2499994605),
Fus3 =>  1/(compartment_1) * (0.0),
Clb5 =>  1/(compartment_1) * (0.0055095267392),
SBF_p =>  1/(compartment_1) * (0.000515873455915),
Swi5_p =>  1/(compartment_1) * (1.10259641365),
Hog1_PP =>  1/(compartment_1) * (0.0),
Clb2_p =>  1/(compartment_1) * (1.13255806709),
Far1_p =>  1/(compartment_1) * (0.0),
Sic1 =>  1/(compartment_1) * (285.604158558),
SBF =>  1/(compartment_1) * (100.0616005253745),
Clb2 =>  1/(compartment_1) * (0.327003199617),
Clb3_Sic1 =>  1/(compartment_1) * (277.885559446),
MBF =>  1/(compartment_1) * (0.0268170651527),
Clb5_Sic1 =>  1/(compartment_1) * (0.89919603461),
Cdc14_p =>  1/(compartment_1) * (338.4307188605),
APC_p =>  1/(compartment_1) * (0.0410087523757),
SBF_Whi5 =>  1/(compartment_1) * (349.937883622),
Whi5_p =>  1/(compartment_1) * (52.41210074),
Cln2 =>  1/(compartment_1) * (200.390588672104),
Whi5 =>  1/(compartment_1) * (620.374730215),
Cdc14 =>  1/(compartment_1) * (3421.56928094),
Far1 =>  1/(compartment_1) * (253.1708654155),
Cln3_Far1_p =>  1/(compartment_1) * (0.0),
Mcm1 =>  1/(compartment_1) * (198.280009914),
Clb2_Sic1 =>  1/(compartment_1) * (22.84051778295),
Clb5_Sic1_Hp =>  1/(compartment_1) * (0.0),
Clb3 =>  1/(compartment_1) * (1.622362084625),
Swi5 =>  1/(compartment_1) * (571.275477765),
Cln2_Far1_p =>  1/(compartment_1) * (0.0),
APC =>  1/(compartment_1) * (47.7610290701),
Swe1 =>  1/(compartment_1) * (16.5140084901),
Swe1_p =>  1/(compartment_1) * (263.70599151),
Mih1 =>  1/(compartment_1) * (107.0),
Sic1_p =>  1/(compartment_1) * (0.612700299515)
]
