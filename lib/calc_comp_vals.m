function comp_vals = calc_comp_vals(params)
lags = [params.tau_2, params.tau_m];
sol = dde23(@(t, n, Z) dde_model_coarse_grained(t, n, Z, params), lags, @(t) dde_history_coarse_grained(t, params),[0, 1000]);
i_last = size(sol.y, 2);

comp_vals = zeros(1, 9);

comp_vals(1) = params.phi_1;
comp_vals(2) = params.delta_1 - params.mu_1 * (1 - params.r_1) * sol.y(3, i_last);
comp_vals(3) = params.mu_1 * params.r_1 * sol.y(3, i_last);
comp_vals(4) = 24 / comp_vals(3);
comp_vals(5) = params.theta_2 * (1 - sol.y(2, i_last) / params.k_2) - (1 - params.r_24) * params.mu_2 * sol.y(3, i_last);
comp_vals(6) = params.r_24 * params.mu_2 * sol.y(3, i_last);
comp_vals(7) = 24 / comp_vals(6);
comp_vals(8) = params.r_24;
comp_vals(9) = params.mu_4 * sol.y(5, i_last); 

end