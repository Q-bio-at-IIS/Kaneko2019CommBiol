%% Parameter estimation of the coarse grained model fixing phi_1 to the value estimated in the DN detail model
% presented in Fig. 4b

addpath(genpath('./'));

data_xray_path = 'data/data_xray.csv';
data_xray = read_data(data_xray_path);

% Parameter names
ps_name = {
    'p_dn' 'omega_dn' 'ndn_0' 'phi_1' 'mu_1' 'delta_1' 'r_1' ...
    'p_dp' 'omega_dp' 'ndp_0'  'theta_2' 'k_2' 'mu_2' 'tau_2' 'r_24' ...
    'p_ctec' 'omega_ctec' 'nctec_0' 'phi_c' 'delta_c' 'mu_c'...
    'p_sp4' 'omega_sp4' 'nsp4_0' 'mu_4' ...
    'p_mtec' 'omega_mtec' 'nmtec_0' 'phi_m' 'phi_m4' 'r_m' 'k_m' 'gamma_mp' 'tau_m' ...
};
cols = {'DN', 'DP', 'cTEC', 'SP4', 'mTEC'};

% Options for parameter estimation
opts=ddeset('Events',@dde_stop_events);
opts.Display = 'iter-detailed';

% Indices of parameters whose bounds are [0, 1]
ratio_ind = [...
    index_of(ps_name, 'p_dn') ...
    index_of(ps_name, 'r_1') ...
    index_of(ps_name, 'p_dp') ...
    index_of(ps_name, 'p_sp4') ...
    index_of(ps_name, 'r_24') ...
    index_of(ps_name, 'p_ctec') ...
    index_of(ps_name, 'p_mtec') ];


% Load initial values of parameters
filename = 'params/coarse_grained_fix_phi.mat';
load(filename)

%% Parameter estimation
% Lower and upper bounds of parameters
lb = ps_give * 0.1;
ub = ps_give * 10;
lb(ratio_ind) = 0;
ub(ratio_ind) = 1.0;

% Fix value of phi_1
est_ind = 1:length(ps_give);
est_ind(est_ind == index_of(ps_name, 'phi_1')) = [];
ps_give_est = ps_give(est_ind);
lb_est = lb(est_ind);
ub_est = ub(est_ind);

[ps_opt_est,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqnonlin(@(ps) diff_model_impulse(subst_est(ps, est_ind, ps_give), data_xray, ps_order, ps_name), ps_give_est, lb_est, ub_est, opts);

ps_opt = ps_give;
ps_opt(est_ind) = ps_opt_est;

save(filename, 'ps_give', 'ps_order', 'ps_opt', 'est_ind', 'lb_est', 'ub_est');

%% 結果のプロットによる確認
params = set_params_order(ps_opt, ps_order, ps_name);
lags = [params.tau_2, params.tau_m];
sol = dde23(@(t, n, Z) dde_model_coarse_grained(t, n, Z, params), lags, @(t) dde_history_coarse_grained(t, params),[0, 49]);
plot_all(sol, cols, data_xray, params);


















