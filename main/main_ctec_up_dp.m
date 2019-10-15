%% Parameter estimation of the possible model (3) Increase of DP by cTEC
% presented in Fig. 5e

addpath(genpath('./'));

data_xray_path = 'data/data_xray.csv';
data_xray = read_data(data_xray_path);

% Parameter names
ps_name = {
    'p_dn' 'omega_dn' 'ndn_0' 'phi_1' 'mu_1' 'delta_1' 'r_1' ...
    'p_dp' 'omega_dp' 'ndp_0'  'r_2' 'k_2' 'mu_2' 'r_24' ...
    'p_ctec' 'omega_ctec' 'nctec_0' 'phi_c' 'delta_c' 'mu_c'...
    'p_sp4' 'omega_sp4' 'nsp4_0' 'mu_4' ...
    'p_mtec' 'omega_mtec' 'nmtec_0' 'phi_m' 'phi_m4' 'r_m' 'k_m' 'gamma_mp' 'tau_m' ...
     };
cols = {'DN', 'DP', 'cTEC', 'SP4', 'mTEC'};

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
load('params/ctec_up_dp.mat');

% Options for parameter estimation
opts=ddeset('Events',@dde_stop_events);
opts.Display = 'iter-detailed';

%% Parameter estimation

% Lower and upper bounds of parameters
lb = ps_give * 0.1;
ub = ps_give * 10;
lb(ratio_ind) = 0;
ub(ratio_ind) = 1.0;

[ps_opt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqnonlin(@(ps) diff_model_handle(@dde_model_ctec_up_dp, ps, data_xray, ps_order, ps_name, {'tau_m'}), ps_give, lb, ub, opts);
save('params/ctec_up_dp.mat', 'ps_give', 'ps_order', 'ps_name', 'ps_opt', 'lb', 'ub');

%% Plot
params = set_params_order(ps_opt, ps_order, ps_name);
lags = [params.tau_m];
sol = dde23(@(t, n, Z) dde_model_ctec_up_dp(t, n, Z, params), lags, @(t) dde_history_coarse_grained(t, params),[0, 49]);
plot_all(sol, cols, data_xray, params);