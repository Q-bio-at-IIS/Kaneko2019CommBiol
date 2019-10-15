function save_boot(ps_boots, boot_datas, boot_filename)
% ps_boots : matrix of bootstrap estimated parameters
% boot_datas : cell array of bootstrap sampling

if exist(boot_filename, 'file')
    load(boot_filename);
    ps_boots_save = cat(1, ps_boots_save, ps_boots);
    boot_datas_save = cat(1, boot_datas_save, boot_datas);
else
    ps_boots_save = ps_boots;
    boot_datas_save = boot_datas;
end
save(boot_filename, 'ps_boots_save', 'boot_datas_save');

end