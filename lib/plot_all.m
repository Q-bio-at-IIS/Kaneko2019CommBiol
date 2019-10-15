function plot_all(sol, cols, data, params)
p_size = 10;

plot_colors.DN = [1.0000e+000   427.4510e-003   713.7255e-003];
plot_colors.cTEC = [7.8431e-003   509.8039e-003   929.4118e-003];
plot_colors.DP = [0.0000e+000     0.0000e+000   666.6667e-003];
plot_colors.SP4 = [541.1765e-003   713.7255e-003    27.4510e-003];
plot_colors.mTEC = [572.5490e-003     0.0000e+000     0.0000e+000];

for i = 1:length(cols)
    col = cols{i};
    n_xray = calc_xray(sol.x, col, params);
    h = semilogy(sol.x,sol.y(i,:)+ n_xray, 'Color', plot_colors.(col));
    hold on;   
    scatter(data.(col).day, data.(col).val, p_size,  get(h,'color'),'filled');
end

xlim([0, 50]);

end