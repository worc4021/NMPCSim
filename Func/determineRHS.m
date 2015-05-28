function [alpha,new] = determineRHS(f,g,old)


tol = 1e-16;

len = length(old)/2;
new = false(size(old));



POSA = [[f(1:len),-ones(len,1)];-1,0;1,0];
POSB = [-g(1:len);0;1];

[POSV,type] = vertexEnumeration(POSA,POSB);
POSV = sortrows(POSV(logical(type),:))';
new(1:len) = POSA(1:len,:)*POSV(:,2)-POSB(1:len)>-tol;
new(old(1:len)) = false;


NEGA = [[f(len+1:end),-ones(len,1)];-1,0;1,0];
NEGB = [-g(len+1:end);0;1];

[NEGV,type] = vertexEnumeration(NEGA,NEGB);
NEGV = sortrows(NEGV(logical(type),:))';
new(len+1:end) = NEGA(1:len,:)*NEGV(:,2)-NEGB(1:len)>-tol;
new(len + find(old(len+1:end)) ) = false;

if and(POSV(1,2)>NEGV(1,2),NEGV(1,2)>tol)
    new(1:len) = old(1:len);
    alpha = NEGV(1,2);
elseif and(POSV(1,2)<NEGV(1,2),POSV(1,2)>tol)
    new(len+1:end) = old(len+1:end);
    alpha = POSV(1,2);
else
    alpha = 1;
end
%     
% 
% 
% 
% ThisStage = [[f,-ones(size(f))];-1,0;1,0];
% thisStage = [-g;0;1];
% [thisV,type] = vertexEnumeration(ThisStage,thisStage);
% idx = logical(type);
% thisV = sortrows(thisV(idx,:))';
% 
% new = ThisStage*thisV(:,2)-thisStage>-tol;
% new = new(1:end-2);
% new(old) = false;
% new(3+len:end) = new(3:len+2);
% if or(thisV(2,2)>thisV(2,1),thisV(1,2)>tol)
%     alpha = thisV(1,2);
% else
%     alpha = thisV(1,end);
% end