function [wVar,wCon,xpVar,xpCon,lambdaVar,lambdaCon,lambdaD,muVar,...
    muCon,muD,zetaVar,zetaCon,zetaD,uVar,uCon,kappaVar,kappaCon,kappaD,...
    rhoVar,rhoCon,rhoD] = forwardsim(wM,wm,xpM,xpm,lambdaM,lambdam,...
    lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,kappaZ,...
    rhoK,rhok,rhoZ,s)

N = s.N;


wVar = zeros(2*N,2);
wCon = zeros(2*N,1);

xpVar = zeros(2*(N+1),2);
xpCon = zeros(2*(N+1),1);

lambdaVar = zeros(N*length(s.g),2);
lambdaCon = zeros(N*length(s.g),1);
lambdaD = zeros(N*length(s.g),1);

muVar = zeros(2*N,2);
muCon = zeros(2*N,1);
muD = zeros(2*N,1);


zetaVar = zeros(1*(N+1),2);
zetaCon = zeros(1*(N+1),1);
zetaD = zeros(1*(N+1),1);

uVar = zeros(N,2);
uCon = zeros(N,1);

kappaVar = zeros(s.nuXi,2);
kappaCon = zeros(s.nuXi,1);
kappaD = zeros(s.nuXi,1);

rhoVar = zeros(N,2);
rhoCon = zeros(N,1);
rhoD = zeros(N,1);

xpVar(1:2,:) = eye(2);

zetaD(1,:) = zetaZ(:,:,1);

curFst = 0;
for i = 1:N
    curLst = curFst + length(s.stage{i}.xi);
    T = xpVar((i-1)*2+1:i*2,:);
    
    uVar(i,:) = uK(:,:,i)*T;
    uCon(i) = uK(:,:,i)*xpCon((i-1)*2+1:i*2) + uk(:,:,i);
    
    kappaVar(curFst+1:curLst,:) = kappaK{i}*T ...
        + kappaZ{i}*zetaVar(i,:);
    kappaCon(curFst+1:curLst) = kappaK{i}*xpCon((i-1)*2+1:i*2) ...
        + kappak{i}+ kappaZ{i}*zetaCon(i);
    kappaD(curFst+1:curLst,:) = kappaZ{i}*zetaD(i,:);
    
    % Finish here!
    
    rhoVar(i,:) = rhoK(:,:,i)*T + rhoZ(:,:,i)*zetaVar(i,:);
    rhoCon(i) = rhoK(:,1:2,i)*xpCon((i-1)*2+1:i*2) ...
        + rhok(:,:,i)+ rhoZ(:,:,i)*zetaCon(i);
    rhoD(i,:) = rhoZ(:,:,i)*zetaD(i,:);
    
    Th = [xpVar((i-1)*2+1:i*2,:);
          uVar(i,:)];
      
    wVar((i-1)*2+1:i*2,:) = wM(:,:,i)*Th;
    wCon((i-1)*2+1:i*2) = wm(:,:,i)+wM(:,:,i)*[xpCon((i-1)*2+1:i*2);uCon(i)];
    
    xpVar(i*2+1:(i+1)*2,:) = xpM(:,:,i)*Th;
    xpCon(i*2+1:(i+1)*2) = xpm(:,:,i)+xpM(:,:,i)*[xpCon((i-1)*2+1:i*2);uCon(i)];
    
    lambdaVar((i-1)*length(s.g)+1:i*length(s.g),:) = lambdaM(:,:,i)*Th ...
        + lambdaZ(:,:,i)*rhoVar(i,:);
    lambdaCon((i-1)*length(s.g)+1:i*length(s.g)) = lambdam(:,:,i) ...
        + lambdaM(:,:,i)*[xpCon((i-1)*2+1:i*2);uCon(i)] ...
        + lambdaZ(:,:,i)*rhoCon(i);
    lambdaD((i-1)*length(s.g)+1:i*length(s.g),:) = lambdaZ(:,:,i)*rhoD(i,:);
    
    muVar((i-1)*2+1:i*2,:) = muM(:,:,i)*Th ...
        + muZ(:,:,i)*rhoVar(i,:);
    muCon((i-1)*2+1:i*2) = mum(:,:,i) ...
        + muM(:,:,i)*[xpCon((i-1)*2+1:i*2);uCon(i)] ...
        + muZ(:,:,i)*rhoCon(i);
    muD((i-1)*2+1:i*2,:) = muZ(:,:,i)*rhoD(i,:);
    
    zetaVar(i,:) = zetaM(:,:,i)*Th ...
        + zetaZ(:,:,i)*rhoVar(i,:);
    zetaCon(i) = zetam(:,:,i) ...
        + zetaM(:,:,i)*[xpCon((i-1)*2+1:i*2);uCon(i)] ...
        + zetaZ(:,:,i)*rhoCon(i);
    zetaD(i,:) = zetaZ(:,:,i)*rhoD(i,:);
    curFst = curLst;
end