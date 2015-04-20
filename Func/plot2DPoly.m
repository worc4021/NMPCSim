function [] = plot2DPoly(A,b)

[VT,typeT] = vertexEnumeration(A,b);
if ~isempty(find(typeT~=1,1))
    VT = VT(logical(typeT),:);
end
order = convhull(VT(:,1),VT(:,2));
fill(VT(order,1),VT(order,2),'r')
