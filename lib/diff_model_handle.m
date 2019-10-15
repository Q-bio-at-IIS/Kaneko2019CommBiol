function diff = diff_model_handle(dde_model, ps_give, data, ps_order, ps_name, lag_names)
params = set_params_order(ps_give, ps_order, ps_name);

d.min = 0;
d.max = 50;

lags = zeros(1, length(lag_names));
for i = 1:length(lag_names)
    lags(i) = params.(lag_names{i});
end

opts = ddeset('Events',@dde_stop_events);
sol = dde23(@(t, n, Z) dde_model(t, n, Z, params), lags, @(t) dde_history_coarse_grained(t, params),[0, 49], opts);

diff = calc_diff_coarse_grained(sol, d, data, params);
end