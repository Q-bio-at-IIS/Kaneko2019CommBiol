function diff = diff_model_impulse(ps_give, data, ps_order, ps_name)
params = set_params_order(ps_give, ps_order, ps_name);

d.min = 0;
d.max = 50;

lags = [params.tau_2, params.tau_m];

opts = ddeset('Events',@dde_stop_events);
sol = dde23(@(t, n, Z) dde_model_coarse_grained(t, n, Z, params), lags, @(t) dde_history_coarse_grained(t, params),[0, 49], opts);

diff = calc_diff_coarse_grained(sol, d, data, params);
end