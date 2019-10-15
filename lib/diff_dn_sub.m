function diff = diff_dn_sub(ps_give, ps_order, data, ps_name)
params = set_params_order(ps_give,ps_order,ps_name);

d.min = 0;
d.max = 50;

lags = [params.tau_2];

sol = dde23(@(t, n, Z) dde_model_dn_sub(t, n, Z, params) ,lags, @(t) dde_history_dn_sub(t, params),[d.min, d.max]);

diff = calc_diff_dn_sub(sol, d, params, data);
end