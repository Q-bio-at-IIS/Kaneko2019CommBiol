function diff = diff_base(sol, d, i, col_name, data, params)
model_day = data.(col_name).day(data.(col_name).day<=d.max & data.(col_name).day>=d.min);
model_vals = deval(sol, model_day);
model_xray = calc_xray(model_day, col_name, params);
model_total = model_vals(i,:).' + model_xray;
diff = log(model_total) - log(data.(col_name).val(data.(col_name).day<=d.max & data.(col_name).day>=d.min));