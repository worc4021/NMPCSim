function [outSys,stepsize] = linesearch(x0,xE,wM,wm,xpM,xpm,lambdaM,lambdam,lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,kappaZ,rhoK,rhok,rhoZ,s)

tol = 1e-9;
xiIdx = logical(genCellConCat(s.stage,'xiIdx'));
wIdx = reshape(s.wIdx,numel(s.wIdx),1);
rhsOld = s.rhsIdx;
rhsMult = double(reshape(rhsOld,numel(rhsOld),1));

[wVar,wCon,xpVar,xpCon,lambdaVar,lambdaCon,lambdaD,muVar,muCon,...
    muD,zetaVar,zetaCon,zetaD,uVar,uCon,kappaVar,kappaCon,kappaD,...
    rhoVar,rhoCon,rhoD] = forwardsim(wM,wm,xpM,xpm,lambdaM,lambdam,...
    lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,...
    kappaZ,rhoK,rhok,rhoZ,s);

wCandNum = kron(eye(s.N),s.G)*(wVar*x0 + wCon)  - diag(rhsMult)*kron(eye(s.N),s.H(:,1:2))*(xpVar(1:end-2,:)*x0+xpCon(1:end-2)) - diag(rhsMult)*kron(eye(s.N),s.H(:,3))*(uVar*x0+uCon) - diag(rhsMult)*repmat(s.g,s.N,1);
wCandDen = (kron(eye(s.N),s.G)*wVar - diag(rhsMult)*kron(eye(s.N),s.H(:,1:2))*xpVar(1:end-2,:) - diag(rhsMult)*kron(eye(s.N),s.H(:,3))*uVar)*(xE-x0);

xiCandNum = - genCellConCat(s.stage,'xi') + cellToBlock(s.stage,'Xi','1:2')*(xpVar(1:end-2,:)*x0+xpCon(1:end-2)) + cellToBlock(s.stage,'Xi','3')*(uVar*x0+uCon);
xiCandDen = (cellToBlock(s.stage,'Xi','1:2')*xpVar(1:end-2,:) + cellToBlock(s.stage,'Xi','3')*uVar)*(xE-x0);

rhsCandNum = kron(eye(s.N),s.H(:,1:2))*(xpVar(1:end-2,:)*x0 + xpCon(1:end-2)) + kron(eye(s.N),s.H(:,3))*(uVar*x0 + uCon) + repmat(s.g,s.N,1);
rhsCandDen = (kron(eye(s.N),s.H(:,1:2))*xpVar(1:end-2,:) + kron(eye(s.N),s.H(:,3))*uVar)*(xE-x0);

rhsStageNum = reshape(rhsCandNum,numel(rhsCandNum)/s.N,s.N);
rhsStageDen = reshape(rhsCandDen,numel(rhsCandDen)/s.N,s.N);

rhsAlpha = zeros(1,s.N);
rhsNew = false(length(s.g),s.N);

for i = 1:s.N
    [rhsAlpha(i),rhsNew(3:end,i)] = determineRHS(rhsStageDen(3:end,i),...
        rhsStageNum(3:end,i),rhsOld(3:end,i));
end

rhsOld(1:2,:) = true(2,s.N);

rhsWholeOld = reshape(rhsOld,numel(rhsOld),1);

distNum = wCandNum - rhsCandNum;
distDen = wCandDen - rhsCandDen;

lambdaCandNum = lambdaVar*x0+lambdaCon;
lambdaCandDen = lambdaVar*(xE-x0);

kappaCandNum = kappaVar*x0+kappaCon;
kappaCandDen = kappaVar*(xE-x0);

wNEQ = and(distDen>tol,rhsWholeOld);
xiNEQ = xiCandDen>tol;
lambdaNEQ = and(lambdaCandDen<-tol,lambdaCandNum>tol);
kappaNEQ = and(kappaCandNum>tol,kappaCandDen<-tol);

xiCandIDX = false(size(xiCandNum));
wCandIDX = false(size(wCandNum));
kappaCandIDX = false(size(kappaCandNum));
lambdaCandIDX = false(size(lambdaCandNum));

[xiCand,temp] = min(-xiCandNum(xiNEQ)./xiCandDen(xiNEQ));
temp2 = find(xiNEQ,temp);
xiCandIDX(temp2(temp)) = true;
[wCand,temp] = min(-distNum(wNEQ)./distDen(wNEQ));
temp2 = find(wNEQ,temp);
wCandIDX(temp2(temp)) = true;
[lambdaCand,temp] = min(-lambdaCandNum(lambdaNEQ)./lambdaCandDen(lambdaNEQ));
if ~isempty(temp)
    temp2 = find(lambdaNEQ,temp);
    lambdaCandIDX(temp2(temp)) = true;
end
[kappaCand,temp] = min(-kappaCandNum(kappaNEQ)./kappaCandDen(kappaNEQ));
if ~isempty(temp)
    temp2 = find(kappaNEQ,temp);
    kappaCandIDX(temp2(temp)) = true;
end

[rhsCand,rhsCandIDX] = min(rhsAlpha);

listed = ones(5,1);
if ~isempty(xiCand)
    listed(1) = xiCand;
end
if ~isempty(wCand)
    listed(2) = wCand;
end
if ~isempty(lambdaCand)
    listed(3) = lambdaCand;
end
if ~isempty(kappaCand)
    listed(4) = kappaCand;
end
if ~isempty(rhsCandIDX)
    listed(5) = rhsCand;
end

[stepsize,sel] = min(listed);

if sel==1
    xiIdx(and(xiNEQ,xiCandIDX)) = ~xiIdx(and(xiNEQ,xiCandIDX));
elseif sel==2
    wIdx(and(wNEQ,wCandIDX)) = ~wIdx(and(wNEQ,wCandIDX));
elseif sel==3
    wIdx(and(lambdaNEQ,lambdaCandIDX)) = ~wIdx(and(lambdaNEQ,lambdaCandIDX));
elseif sel==4
    xiIdx(and(kappaNEQ,kappaCandIDX)) = ~xiIdx(and(kappaNEQ,kappaCandIDX));
elseif sel==5
    rhsOld(3:end,rhsCandIDX) = rhsNew(3:end,rhsCandIDX);
else
    stepsize = 1e10;
end

outSys = s;

for i = 1:s.N
    outSys.stage{i}.xiIdx = xiIdx(1:length(s.stage{i}.xi));
    xiIdx(1:length(s.stage{i}.xi)) = [];
end
outSys.wIdx = reshape(wIdx,numel(wIdx)/s.N,s.N);
outSys.rhsIdx = rhsOld;