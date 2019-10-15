%% Bootstrap estimation of the coarse grained model
% presented in Fig. 2c, Fig. 3, Supplementary Fig. 2, and Supplementary Fig. 3

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

load('params/coarse_grained.mat');

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


%% Calculate variance
params = set_params_order(ps_opt, ps_order, ps_name);
lags = [params.tau_2, params.tau_m];
sol = dde23(@(t, n, Z) dde_model_coarse_grained(t, n, Z, params), lags, @(t) dde_history_coarse_grained(t, params),[0, 49]);
d.min = 0;
d.max = 50;

for i = 1:length(cols)
    col = cols{i};
    diff_est.(col) = diff_base(sol, d, i, col, data_xray, params);
    sigma_est.(col) = sqrt(diff_est.(col).' * diff_est.(col) / (length(diff_est.(col))-1));
end

%% Bootstrap estimation

% Lower and upper bounds of parameters
lb = ps_opt * 0.1;
ub = ps_opt * 10;
lb(ratio_ind) = 0;
ub(ratio_ind) = 1.0;

boot_filename = './params/bootstrap.mat';
n_iter = 10;

for j = 1:100
    disp(['%%%%%%%%%%%%%%%%%% ' num2str(j) ' %%%%%%%%%%%%%%%%%%']);
    rng('shuffle')

    ps_boots = zeros(n_iter, length(ps_opt));
    boot_datas = cell(n_iter, 1);

    tic
    for i = 1:n_iter
        disp([num2str(i) '/' num2str(n_iter) ]);    

        % Generate bootstrap sample
        boot_datas{i,1} = boot_sample(sol, params, sigma_est, data_xray, cols);
        
        % Bootstrap estimation
        [ps_boots(i,:), resnorm, residual, exitflag, output, lambda, jacobian] = ...
            lsqnonlin(@(ps) diff_model_impulse(ps, boot_datas{i,1}, ps_order, ps_name), ps_opt, lb, ub, opts);
    end
    toc
    save_boot(ps_boots, boot_datas, boot_filename); 
end

%% Plot
load('params/bootstrap.mat');

for i = 1:100
    params = set_params_order(ps_boots_save(i,:), ps_order, ps_name);
    lags = [params.tau_2, params.tau_m];
    sol = dde23(@(t, n, Z) dde_model_coarse_grained(t, n, Z, params), lags, @(t) dde_history_coarse_grained(t, params),[0, 49]);
    plot_all(sol, cols, data_xray, params);
end
