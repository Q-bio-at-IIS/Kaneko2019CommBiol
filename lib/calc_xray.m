function val_xray = calc_xray(days, col, params)
if strcmp(col,'DN')
    val_xray =  params.ndn_0 * (1-params.p_dn) * exp(-params.omega_dn * days);
elseif strcmp(col,'DP')
    val_xray =  params.ndp_0 * (1-params.p_dp) * exp(-params.omega_dp * days);   
elseif strcmp(col,'SP4')
    val_xray =  params.nsp4_0 * (1-params.p_sp4) * exp(-params.omega_sp4 * days);
elseif strcmp(col,'cTEC')
    val_xray =  params.nctec_0 * (1-params.p_ctec) * exp(-params.omega_ctec * days);
elseif strcmp(col,'mTEC')
    val_xray =  params.nmtec_0 * (1-params.p_mtec) * exp(-params.omega_mtec * days);
end

end