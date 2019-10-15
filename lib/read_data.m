function data = read_data(data_path)
df = readtable(data_path);
data = struct();
day = df.Properties.VariableNames{1};
for i = 1:length(df.Properties.VariableNames)
    col = df.Properties.VariableNames{i};
    data.(col) = table(df.(day)(~isnan(df.(col))),  df.(col)(~isnan(df.(col))), 'VariableNames', {'day','val'});
end

end