function dndt = dde_model_no_mtec_self(t, n, Z, params)
DN = n(1);
DP = n(2);
cTEC = n(3);
SP4 = n(4);
mTEC = n(5);
cTEC_tau = Z(3,1);
DP_tau = Z(2,2);

next_DN = params.phi_1 - params.mu_1 * cTEC .* DN + params.delta_1 * DN;
next_DP = params.mu_1*params.r_1*cTEC*DN + params.theta_2 *DP*(1-DP/params.k_2) - params.mu_2 * cTEC_tau * DP;
next_cTEC =  params.phi_c - (params.delta_c - params.mu_c * DN) .* cTEC;
next_SP4 =  params.mu_2 * params.r_24 * DP * cTEC_tau - params.mu_4 * mTEC * SP4;
next_mTEC = params.phi_m  + params.phi_m4 * SP4 + (-params.gamma_mp * DP_tau - params.mu_m) * mTEC;

dndt = [next_DN; next_DP; next_cTEC; next_SP4; next_mTEC];