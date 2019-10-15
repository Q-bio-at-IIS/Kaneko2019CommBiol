function ind = index_of(cell, str)
indCont = strcmp(cell, str);
ind = find(indCont == 1);
end