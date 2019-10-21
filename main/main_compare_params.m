addpath(genpath('./'));

%% initialize
rows= {...
    'DN influx', 'DN apparent growth', 'DN diff', 'DN stay', ...
    'DP apparent growth',  'DP diff to SP4', 'DP stay', 'DP diff to SP4(%)'...
    'SP4 apparent outflux', ...    
    };
n_row = length(rows);
table_params = table(zeros(n_row,1), zeros(n_row,1), zeros(n_row,1),zeros(n_row,1),zeros(n_row,1),zeros(n_row,1),zeros(n_row,1),zeros(n_row,1), zeros(n_row,1),  ...
    'RowNames',rows, ...
    'VariableNames', {'ThomasVaslin2008','Manesso2012','Sinclair2013','Sawicka2014' 'Moleriu2014',  'Zaharie2016', 'Ours_opt', 'Ours_lb','Ours_ub' });

rates = struct();



%% ThomasVaslin2008
ps_thomas = struct();

ps_thomas.delta = 0;
ps_thomas.sigma_N = 0.02 * 10^6;
ps_thomas.n = 127;
% DN
ps_thomas.p_N = 0.23;
ps_thomas.Ns = [0.09*10^6 0.18*10^6 0.36*10^6 0.72*10^6 3*10^3];
ps_thomas.n_N = 4;
ps_thomas.alpha_muN = 0.29;
ps_thomas.mu_NN_4 = 0.32 * 10^6;

% DP
ps_thomas.Ps = [0.07 0.14 0.29 0.58 1.16 0.94 0.08 40.29] * 10^6;
ps_thomas.n_P = 7;
ps_thomas.p_P = 4.5;
ps_thomas.alpha_4 = 0.06;
ps_thomas.alpha_8 = 0.01;
ps_thomas.mu_LP = 0.37;

% SP
ps_thomas.S4s = [4.02 2.49 0.01] * 10^6;
ps_thomas.p_S = 0.23;
ps_thomas.alpha_e = 0.99;
ps_thomas.eS4_1= 1.26 * 10^6;
ps_thomas.eS4_2 = 1.13 * 10^6;

% rate
rates.thomas = struct();
% DN
rates.thomas.DN_inflow = ps_thomas.sigma_N;
rates.thomas.DN_growth = ps_thomas.p_N * sum(ps_thomas.Ns(1:4)) / sum(ps_thomas.Ns);
rates.thomas.DN_death = ps_thomas.delta;
rates.thomas.DN_diff = ps_thomas.mu_NN_4 / sum(ps_thomas.Ns);
rates.thomas.DN_ts = zeros(1,5);
rates.thomas.DN_ts(1) = 1 / (ps_thomas.p_N + ps_thomas.delta);
for i = 2:5
    rates.thomas.DN_ts(i) = ps_thomas.p_N * rates.thomas.DN_ts(i-1)  / (ps_thomas.p_N + ps_thomas.delta + min(( ps_thomas.alpha_muN * i-1)^ps_thomas.n,100) );
end
rates.thomas.DN_stay = 17.6*24;

% DP
rates.thomas.DP_growth = ps_thomas.p_P * sum(ps_thomas.Ps(1:7)) / sum(ps_thomas.Ps);
rates.thomas.DP_death = (ps_thomas.mu_LP * (1 - ps_thomas.alpha_4 - ps_thomas.alpha_8) * ps_thomas.Ps(ps_thomas.n_P+1)) / sum(ps_thomas.Ps);
rates.thomas.DP_diff4 = (ps_thomas.mu_LP * ps_thomas.alpha_4 * ps_thomas.Ps(ps_thomas.n_P+1)) /  sum(ps_thomas.Ps);
rates.thomas.DP_diff8 = (ps_thomas.mu_LP * ps_thomas.alpha_8* ps_thomas.Ps(ps_thomas.n_P+1)) / sum(ps_thomas.Ps);
rates.thomas.DP_stay = (1.2+2.7)*24;
rates.thomas.DP_to_4_ratio = ps_thomas.alpha_4;

% SP4
rates.thomas.SP4_growth = ps_thomas.p_S * sum(ps_thomas.S4s(1:2)) / sum(ps_thomas.Ps);
rates.thomas.SP4_death = ps_thomas.delta;
rates.thomas.SP4_diff = (ps_thomas.eS4_1 + ps_thomas.eS4_2) / sum(ps_thomas.Ps);
rates.thomas.SP4_stay = 5.8*24;

% table
table_params('DN influx', 'ThomasVaslin2008') = {rates.thomas.DN_inflow};
table_params('DN apparent growth', 'ThomasVaslin2008') = {rates.thomas.DN_growth - rates.thomas.DN_death};
table_params('DN diff', 'ThomasVaslin2008') =  {rates.thomas.DN_diff};
table_params('DN stay', 'ThomasVaslin2008') = {rates.thomas.DN_stay};
table_params('DP apparent growth', 'ThomasVaslin2008') = {rates.thomas.DP_growth-rates.thomas.DP_death-rates.thomas.DP_diff8};
table_params('DP diff to SP4', 'ThomasVaslin2008') = {rates.thomas.DP_diff4};
table_params('DP stay', 'ThomasVaslin2008') = {rates.thomas.DP_stay};
table_params('DP diff to SP4(%)', 'ThomasVaslin2008') = {rates.thomas.DP_to_4_ratio};
table_params('SP4 apparent outflux', 'ThomasVaslin2008') = {-rates.thomas.SP4_growth+rates.thomas.SP4_death+rates.thomas.SP4_diff};

%% Sawicka2014
ps_sawicka = struct();
ps_sawicka.n_1 = 88.39 * 10^6;
ps_sawicka.n_2 = 8.45 * 10^6;
ps_sawicka.mu_1  = 0.263;
ps_sawicka.mu_2 = 1.372;
ps_sawicka.varphi_4 = 0.140;
ps_sawicka.varphi_8 = 0.134;
ps_sawicka.tau_1 = 60;
ps_sawcika.tau_2 = 16;
ps_sawicka.n_4 = 12.33 * 10^6;
ps_sawicka.tau_4 = 96;
ps_sawicka.lambda_4 = 0.181;
ps_sawicka.mu_4 = 0.04;
ps_sawicka.xi_4 = 0.231;

% DP
rates.sawicka = struct();
rates.sawicka.DP_death = (ps_sawicka.mu_1 * ps_sawicka.n_1 + ps_sawicka.mu_2 * ps_sawicka.n_2) / (ps_sawicka.n_1 + ps_sawicka.n_2);
rates.sawicka.DP_diff4 = ps_sawicka.varphi_4 * ps_sawicka.n_2 / (ps_sawicka.n_1 + ps_sawicka.n_2);
rates.sawicka.DP_diff8 = ps_sawicka.varphi_8 * ps_sawicka.n_2 / (ps_sawicka.n_1 + ps_sawicka.n_2);
rates.sawicka.DP_stay = ps_sawicka.tau_1 + ps_sawcika.tau_2;
rates.sawicka.DP_to_4_ratio = rates.sawicka.DP_diff4 / (rates.sawicka.DP_diff4+rates.sawicka.DP_diff8+rates.sawicka.DP_stay);

% SP4
rates.sawicka.SP_growth = ps_sawicka.lambda_4;
rates.sawicka.SP_death = ps_sawicka.mu_4;
rates.sawicka.SP_diff = ps_sawicka.xi_4;
rates.sawicka.SP_stay = ps_sawicka.tau_4;

% table
table_params('DP apparent growth', 'Sawicka2014') = {-rates.sawicka.DP_death - rates.sawicka.DP_diff8};
table_params('DP diff to SP4', 'Sawicka2014')  = {rates.sawicka.DP_diff4};
table_params('DP stay', 'Sawicka2014') = {rates.sawicka.DP_stay};
table_params('DP diff to SP4(%)', 'Sawicka2014') = {rates.sawicka.DP_to_4_ratio};
table_params('SP4 apparent outflux', 'Sawicka2014') = {-rates.sawicka.SP_growth+rates.sawicka.SP_diff+rates.sawicka.SP_death};

%% Moleriu
ps_moleriu = struct();
ps_moleriu.N = 4.5 * 10^7;
ps_moleriu.P = 10.1 * 10^7;
ps_moleriu.M_4 = 2.0 * 10^7;
ps_moleriu.M_8 = 1.4 * 10.^7;
ps_moleriu.Z = ps_moleriu.N+ps_moleriu.P+ps_moleriu.M_4+ps_moleriu.M_8;
ps_moleriu.K  = 32*10^7;
ps_moleriu.s_N = 0.45;
ps_moleriu.r_P = 0.41;
ps_moleriu.d_P = 0.30;
ps_moleriu.s_4 = 0.03;
ps_moleriu.s_8 = 0.04;
ps_moleriu.r_4 = 0.31;
ps_moleriu.d_4_s_o4 = 0.28;

% rate
rates.moleriu = struct();

% DN
rates.moleriu.DN_diff = ps_moleriu.s_N;

% DP
rates.moleriu.DP_growth_death = ps_moleriu.r_P * (1- ps_moleriu.Z / ps_moleriu.K) - ps_moleriu.d_P;
rates.moleriu.DP_diff4 = ps_moleriu.s_4;
rates.moleriu.DP_diff8 = ps_moleriu.s_8;
rates.moleriu.DP_to_4_ratio =  ps_moleriu.s_4 / (ps_moleriu.s_4 +  ps_moleriu.s_8 + ps_moleriu.d_P);

% SP4
rates.moleriu.SP4_growth_death_diff = -ps_moleriu.r_4 * (1 - ps_moleriu.Z / ps_moleriu.K ) + ps_moleriu.d_4_s_o4;

% table
table_params('DN diff', 'Moleriu2014') = {rates.moleriu.DN_diff};
table_params('DP apparent growth', 'Moleriu2014') = {rates.moleriu.DP_growth_death - rates.moleriu.DP_diff8};
table_params('DP diff to SP4', 'Moleriu2014') = {rates.moleriu.DP_diff4};
table_params('DP diff to SP4(%)', 'Moleriu2014') = {rates.moleriu.DP_to_4_ratio};
table_params('SP4 apparent outflux', 'Moleriu2014') = {rates.moleriu.SP4_growth_death_diff};



%% Zaharie
ps_zaharie = struct();
ps_zaharie.t_0 = mean([9.99, 10.05]);
ps_zaharie.b_0  = mean([0.0012859, 0.00131817]) * 10^7;
ps_zaharie.beta = mean([2.24, 2.28]);
ps_zaharie.tau_b = mean([14.19, 14.24]);
ps_zaharie.s_N = mean([0.027, 0.029]);
ps_zaharie.b_N = mean([0.06, 0.07]);
ps_zaharie.c_N = mean([5.88, 5.98]);
ps_zaharie.d_N = mean([0.0003, 0.001]);
ps_zaharie.b_P = mean([0.08, 0.09]);
ps_zaharie.c_P = mean([4.26, 4.34]);
ps_zaharie.d_P = mean([0.0021, 0.0034]);
ps_zaharie.s_4 = mean([0.017, 0.018]);
ps_zaharie.s_8  = mean([0.006, 0.007]);
ps_zaharie.b_4 = mean([0.11, 0.12]);
ps_zaharie.c_48 = mean([6.18, 6.29]);
ps_zaharie.d_4 = mean([0.16, 0.17]);
ps_zaharie.s_o4 =  mean([0.003, 0.006]);
t = 7*7 + 49;

% rates
rates.zaharie = struct();

% DN
rates.zaharie.DN_inflow = ps_zaharie.b_0 / (1 + exp(- ps_zaharie.beta * ( t - ps_zaharie.tau_b )) );
rates.zaharie.DN_growth_death = ps_zaharie.b_N * ps_zaharie.c_N * exp(-ps_zaharie.b_N * (t - ps_zaharie.t_0)) - ps_zaharie.d_N;
rates.zaharie.DN_diff = ps_zaharie.s_N;

% DP
rates.zaharie.DP_growth_death =  ps_zaharie.b_P * ps_zaharie.c_P * exp(-ps_zaharie.b_P * (t - ps_zaharie.t_0)) - ps_zaharie.d_P;
rates.zaharie.DP_diff4 = ps_zaharie.s_4;
rates.zaharie.DP_diff8 = ps_zaharie.s_8;
rates.zaharie.DP_to_4_ratio = ps_zaharie.s_4 / (ps_zaharie.s_4+ps_zaharie.s_8+ps_zaharie.d_P);

% SP
rates.zaharie.SP4_growth_death_diff = -ps_zaharie.b_4 * ps_zaharie.c_48 * exp(- ps_zaharie.b_4 * (t-ps_zaharie.t_0) ) + ps_zaharie.d_4 + ps_zaharie.s_o4;

% table
table_params('DN influx' , 'Zaharie2016') =  {rates.zaharie.DN_inflow};
table_params('DN apparent growth', 'Zaharie2016')  = {rates.zaharie.DN_growth_death};
table_params('DN diff', 'Zaharie2016') = {rates.zaharie.DN_diff};
table_params('DP apparent growth', 'Zaharie2016') = {rates.zaharie.DP_growth_death - rates.zaharie.DP_diff8};
table_params('DP diff to SP4', 'Zaharie2016') = {rates.zaharie.DP_diff4};
table_params('DP diff to SP4(%)', 'Zaharie2016') = {rates.zaharie.DP_to_4_ratio};
table_params('SP4 apparent outflux', 'Zaharie2016') = {rates.zaharie.SP4_growth_death_diff};


%% Manesso
ps_manesso = struct();
ps_manesso.N_DN1pre = 1.13 * 10^2;
ps_manesso.N_DN1 = 2.73 * 10^4;
ps_manesso.N_DN1_11 = 1.23 * 10^4;
ps_manesso.N_DN2 = 2.88 * 10^4;
ps_manesso.N_DN3 = 2.63 * 10^6;
ps_manesso.N_pDP = 1.57 * 10^6;
ps_manesso.tau = 2.27;
ps_manesso.c_DN1_11 = 1.0;
ps_manesso.T_DN1 = 1.10;
ps_manesso.T_DN2 = 0.631;
ps_manesso.T_DN3 = 2.65;
ps_manesso.T_pDP = 0.514;
ps_manesso.c_DN2 = 0.58;
ps_manesso.c_DN3 = 0.216;
ps_manesso.c_pDP = 0.480;
ps_manesso.d_DN2 = 0.028;
ps_manesso.d_DN3 = 0.216;
ps_manesso.d_pDP = 0.073;
ps_manesso.MTT_DN1 = 9.76*24;
ps_manesso.MTT_DN2 = 2.47*24;
ps_manesso.MTT_DN3 = 1.68*24;
ps_manesso.MTT_pDP = 0.599*24;
ps_manesso.N = ps_manesso.N_DN1+ps_manesso.N_DN2+ps_manesso.N_DN3+ps_manesso.N_pDP;

% rate
rates.manesso = struct();
rates.manesso.DN_inflow = ps_manesso.N_DN1pre / ps_manesso.tau;
rates.manesso.DN_growth = ((ps_manesso.N_DN1 - ps_manesso.N_DN1_11) / ps_manesso.T_DN1 + ...
    (1-ps_manesso.c_DN2-ps_manesso.d_DN2) * ps_manesso.N_DN2 / ps_manesso.T_DN2 + ...
    (1-ps_manesso.c_DN3-ps_manesso.d_DN3) * ps_manesso.N_DN3 / ps_manesso.T_DN3 + ...    
    (1-ps_manesso.c_pDP-ps_manesso.d_pDP) * ps_manesso.N_pDP / ps_manesso.T_pDP ...
    ) / ps_manesso.N;

rates.manesso.DN_death = ( ...
    ps_manesso.d_DN2 * ps_manesso.N_DN2 / ps_manesso.T_DN2 + ...
    ps_manesso.d_DN3 * ps_manesso.N_DN3 / ps_manesso.T_DN3 + ...    
    ps_manesso.d_pDP * ps_manesso.N_pDP / ps_manesso.T_pDP ...
    ) / ps_manesso.N;

rates.manesso.DN_diff = ps_manesso.c_pDP * ps_manesso.N_pDP / (ps_manesso.T_pDP * ps_manesso.N);
rates.manesso.DN_stay = ps_manesso.MTT_DN1 + ps_manesso.MTT_DN2 + ps_manesso.MTT_DN3 + ps_manesso.MTT_pDP;

% table
table_params('DN influx' , 'Manesso2012') =  {rates.manesso.DN_inflow};
table_params('DN apparent growth', 'Manesso2012')  = {rates.manesso.DN_growth - rates.manesso.DN_death};
table_params('DN diff', 'Manesso2012') = {rates.manesso.DN_diff};
table_params('DN stay', 'Manesso2012') = {rates.manesso.DN_stay};


%% sinclair
table_params('DP stay', 'Sinclair2013') = {5 * 24};

%% Our
load('params/coarse_grained.mat')
load('params/bootstrap.mat')

ps_name = {
    'p_dn' 'omega_dn' 'ndn_0' 'phi_1' 'mu_1' 'delta_1' 'r_1' ...
    'p_dp' 'omega_dp' 'ndp_0'  'theta_2' 'k_2' 'mu_2' 'tau_2' 'r_24' ...
    'p_ctec' 'omega_ctec' 'nctec_0' 'phi_c' 'delta_c' 'mu_c'...
    'p_sp4' 'omega_sp4' 'nsp4_0' 'mu_4' ...
    'p_mtec' 'omega_mtec' 'nmtec_0' 'phi_m' 'phi_m4' 'r_m' 'k_m' 'gamma_mp' 'tau_m' ...
};

params = set_params_order(ps_opt, ps_order, ps_name);
comp_vals_opt = calc_comp_vals(params);

n_comp_val = 9;
alpha = 0.05;

ci_params_prcs = zeros(n_comp_val,2);
for i = 1:n_comp_val
    comp_val_sorted = sort(comp_vals{:,i});
    ci_params_prcs(i,1) = prctile(comp_val_sorted, 100*alpha/2);    
    ci_params_prcs(i,2) = prctile(comp_val_sorted, 100*(1-alpha/2));
end

% table
table_params('DN influx', 'Ours_opt') = {comp_vals_opt(1)};
table_params('DN apparent growth', 'Ours_opt') = {comp_vals_opt(2)};
table_params('DN diff', 'Ours_opt') =  {comp_vals_opt(3)};
table_params('DN stay', 'Ours_opt') = {comp_vals_opt(4)};
table_params('DP apparent growth', 'Ours_opt') = {comp_vals_opt(5)};
table_params('DP diff to SP4', 'Ours_opt') = {comp_vals_opt(6)};
table_params('DP stay', 'Ours_opt') = {comp_vals_opt(7)};
table_params('DP diff to SP4(%)', 'Ours_opt') = {comp_vals_opt(8)};
table_params('SP4 apparent outflux', 'Ours_opt') = {comp_vals_opt(9)};

table_params('DN influx', 'Ours_lb') = {ci_params_prcs(1,1)};
table_params('DN apparent growth', 'Ours_lb') = {ci_params_prcs(2,1)};
table_params('DN diff', 'Ours_lb') =  {ci_params_prcs(3,1)};
table_params('DN stay', 'Ours_lb') = {ci_params_prcs(4,1)};
table_params('DP apparent growth', 'Ours_lb') = {ci_params_prcs(5,1)};
table_params('DP diff to SP4', 'Ours_lb') = {ci_params_prcs(6,1)};
table_params('DP stay', 'Ours_lb') = {ci_params_prcs(7,1)};
table_params('DP diff to SP4(%)', 'Ours_lb') = {ci_params_prcs(8,1)};
table_params('SP4 apparent outflux', 'Ours_lb') = {ci_params_prcs(9,1)};

table_params('DN influx', 'Ours_ub') = {ci_params_prcs(1,2)};
table_params('DN apparent growth', 'Ours_ub') = {ci_params_prcs(2,2)};
table_params('DN diff', 'Ours_ub') =  {ci_params_prcs(3,2)};
table_params('DN stay', 'Ours_ub') = {ci_params_prcs(4,2)};
table_params('DP apparent growth', 'Ours_ub') = {ci_params_prcs(5,2)};
table_params('DP diff to SP4', 'Ours_ub') = {ci_params_prcs(6,2)};
table_params('DP stay', 'Ours_ub') = {ci_params_prcs(7,2)};
table_params('DP diff to SP4(%)', 'Ours_ub') = {ci_params_prcs(8,2)};
table_params('SP4 apparent outflux', 'Ours_ub') = {ci_params_prcs(9,2)};
