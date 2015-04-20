function out = genCellConCat(inp,character)

out = [];
for i = 1:length(inp)
    out = [out;eval(['inp{i}.',character])];
end