function params = set_params_order(ps_val, ps_order, ps_name_give)
global ps_name;

if nargin == 2
    ps_name_use = ps_name;
else
    ps_name_use = ps_name_give;
end

if length(ps_val) ~= length(ps_name_use)
    error('parameter name length not match');
end
params = struct();
for i = 1:length(ps_val)  
    params.(ps_name_use{i}) = ps_val(i) * ps_order(i);
end
end