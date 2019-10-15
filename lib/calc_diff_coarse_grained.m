function diff = calc_diff_coarse_grained(sol, d, data, params)
if sol.x(end) < 49
    len_data = length(data.DN.day(data.DN.day<=d.max & data.DN.day>=d.min))...
    +length(data.cTEC.day(data.cTEC.day<=d.max & data.cTEC.day>=d.min))...
    +length(data.DP.day(data.DP.day<=d.max & data.DP.day>=d.min))...
    +length(data.SP4.day(data.SP4.day<=d.max & data.SP4.day>=d.min))...
    +length(data.mTEC.day(data.mTEC.day<=d.max & data.mTEC.day>=d.min));
    diff(1:len_data) = 1.0E+8+rand(1)*1000000;
else
    cols = char('DN', 'DP', 'cTEC', 'SP4', 'mTEC');
    diff = [];
    for i = 1:size(cols, 1)
        col = strtrim(cols(i,:));
        diff = cat(1, diff, diff_base(sol, d, i, col, data, params));
    end
end
end