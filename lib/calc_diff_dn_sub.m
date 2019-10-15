function diff = calc_diff_dn_sub(sol, d, params, df_notnull)

if sol.x(end) < d.max - 0.5

    len_df_notnull = length(df_notnull.DN.day(df_notnull.DN.day<=d.max & df_notnull.DN.day>=d.min)) * 5 ...
    +length(df_notnull.cTEC.day(df_notnull.cTEC.day<=d.max & df_notnull.cTEC.day>=d.min));
    
    diff(1:len_df_notnull) = 1.0E+8;
else

    % xray params
    n0xs = [params.ndn10 * (1-params.pdn1); ...
        params.ndn20 * (1-params.pdn2); ...
        params.ndn30 * (1-params.pdn3); ...
        params.ndn40 * (1-params.pdn4); ...
        params.nctec0 * (1-params.pctec); ...
        params.ndp0 * (1-params.pdp) ...
    ];
    omegas = [params.omegadn1 params.omegadn2 params.omegadn3 params.omegadn4 params.omegactec params.omegadp];
    
    % DN
    model_dn_day = df_notnull.DN.day(df_notnull.DN.day<=d.max & df_notnull.DN.day>=d.min);
    model_dn_vals = deval(sol, model_dn_day);
    xrays_dn = n0xs .* exp(-omegas' *  model_dn_day');
    ntots_dn = model_dn_vals + xrays_dn;
    

    diff_dn = zeros(1, 4 * length(model_dn_day));    
    cols_dnsub = {'DN1', 'DN2', 'DN3', 'DN4'};
    for i = 1:4
        col = cols_dnsub{i};
        diff_dn( (i-1) * length(model_dn_day) + 1: i * length(model_dn_day) ) = log(ntots_dn(i,:)) - log(df_notnull.(col).val(df_notnull.(col).day<=d.max & df_notnull.(col).day>=d.min)).';
    end    


    % cTEC
    model_ctec_day = df_notnull.cTEC.day(df_notnull.cTEC.day<=d.max & df_notnull.cTEC.day>=d.min);
    model_ctec_vals = deval(sol, model_ctec_day);
    xrays_ctec = n0xs .* exp(-omegas' *  model_ctec_day');
    ntots_ctec = model_ctec_vals + xrays_ctec;    
    diff_ctec = log(ntots_ctec(5,:)) - log(df_notnull.cTEC.val(df_notnull.cTEC.day<=d.max & df_notnull.cTEC.day>=d.min)).';
    
    % DP
    model_dp_day = df_notnull.DP.day(df_notnull.DP.day<=d.max & df_notnull.DP.day>=d.min);
    model_dp_vals = deval(sol, model_dp_day);
    xrays_dp = n0xs .* exp(-omegas' *  model_dp_day');
    ntots_dp = model_dp_vals + xrays_dp;    
    diff_dp = log(ntots_dp(6,:)) - log(df_notnull.DP.val(df_notnull.DP.day<=d.max & df_notnull.DP.day>=d.min)).';
    
    diff = cat(2, diff_dn, diff_ctec, diff_dp).';
end
end