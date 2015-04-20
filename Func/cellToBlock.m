function out = cellToBlock(inp,character,varargin)

out = [];

if length(varargin)<1
    range = ':';
else
    range = varargin{1};
end

for i = 1:length(inp)
    out = blkdiag(out,eval(['inp{i}.',character,'(:,',range,')']));
end