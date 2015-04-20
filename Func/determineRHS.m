function [alpha,new] = determineRHS(f,g,old)


tol = 1e-16;

ThisStage = [[f,-ones(size(f))];-1,0;1,0];
thisStage = [-g;0;1];
[thisV,type] = vertexEnumeration(ThisStage,thisStage);
idx = logical(type);
thisV = sortrows(thisV(idx,:))';

new = ThisStage*thisV(:,2)-thisStage>-tol;
new = new(1:end-2);
new(old) = false;
if or(thisV(2,2)>thisV(2,1),thisV(1,2)>tol)
    alpha = thisV(1,2);
else
    alpha = thisV(1,end);
end