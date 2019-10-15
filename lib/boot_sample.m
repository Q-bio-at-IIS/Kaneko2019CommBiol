function boot_data = boot_sample(sol, params, sigma_est, data, cols)
% simulation using 'ps_est'
boot_data = {};
for i = 1:length(cols)
    col = cols{i};
    % Get sample of same size with original data
    boot_day = data.(col).day;
    boot_vals_base = deval(sol, boot_day);
    boot_vals_xray = calc_xray(boot_day, col, params);
    % Add noise of variation 'sigma_est'
    boot_val = (boot_vals_base(i,:) + boot_vals_xray.') .* exp(randn(1,length(boot_vals_base(i,:)))*sigma_est.(col));        
    boot_data.(col) = table(boot_day, boot_val.', 'VariableNames',{'day','val'});
end
end