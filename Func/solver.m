function [wM,wm,xpM,xpm,lambdaM,lambdam,lambdaZ,muM,mum,muZ,zetaM,zetam,...
    zetaZ,uK,uk,kappaK,kappak,kappaZ,rhoK,rhok,rhoZ] = solver(s)

actW = s.wIdx;

qc = zeros(2,1);
Pc = s.P;
Ch = zeros(0,3);
W = blkdiag(s.Q,s.R);
phih = [];

initMats;

for n = s.N:-1:1
    s.m = n;
    nlambda = numel(find(actW(:,n)));
    sol = maxWrapper(Pc,qc,Ch,phih,s);
    
    wM(:,:,n) = sol.wM;
    wm(:,1,n) = sol.wm;
    
    xpM(:,:,n) = sol.xpM;
    xpm(:,1,n) = sol.xpm;
    
    defect = size(sol.lambdaZ,2);
    
    lambdaM(actW(:,n),:,n) = sol.lambdaM;
	lambdam(actW(:,n),1,n) = sol.lambdam;
    lambdaZ(actW(:,n),1:defect,n) = sol.lambdaZ;

    muM(:,:,n) = sol.muM;
    mum(:,1,n) = sol.mum;
    muZ(:,1:defect,n) = sol.muZ;
    
    nzeta = length(sol.zetam);
    zetaM(1:nzeta,:,n) = sol.zetaM;
    zetam(1:nzeta,1,n) = sol.zetam;
    zetaZ(1:nzeta,1:defect,n) = sol.zetaZ;
    
    T = xpM(:,:,n);
	mh = xpm(:,1,n);
    
    Ph = T'*Pc*T-s.gammasq*wM(:,:,n)'*wM(:,:,n);
    qh = (qc'*T+mh'*Pc*T-s.gammasq*wm(:,1,n)'*wM(:,:,n))';
    
    Z = [sol.lambdaZ;
         sol.muZ;
         sol.zetaZ];
    
    if ~isempty(Z)
        C = -Z'*[s.A,s.B;
                 s.stage{n}.Xi(actW(n,:),:);
                 zeros(size(Ch,1),s.nX+s.nU),-Ch(s.nX+s.nU:end)];
        phi = Z'*[zeros(nlambda,1);
                  s.g(actW(n,:));
                  phih];
    else
        C = zeros(0,3);
        phi = [];
    end
    
    actX = s.stage{n}.xiIdx;
    
	nkappa = numel(find(actX));
    
    if nkappa>0
        0;
    end
    
    sol = minWrapper(Ph,qh,C,phi,actX,s);
    
    uK(:,:,n) = sol.uK;
    uk(:,:,n) = sol.uk;
    
    defect = max(size(sol.kappaZ,2),size(sol.rhoZ,2));
    
    
    kappaK{n} = zeros(length(s.stage{n}.xi),2);
    kappaK{n}(actX,:) = sol.kappaK;
    kappak{n} = zeros(length(s.stage{n}.xi),1);
    kappak{n}(actX,:) = sol.kappak;
    kappaZ{n} = zeros(length(s.stage{n}.xi),1);
    kappaZ{n}(actX,1:defect) = sol.kappaZ;
    
    nrho = size(sol.rhoK,1);
    rhoK(1:nrho,:,n) = sol.rhoK;
    rhok(1:nrho,:,n) = sol.rhok;
    rhoZ(1:nrho,1:defect,:) = sol.rhoZ;
    
    Z = [sol.kappaZ;
         sol.rhoZ];
	if ~isempty(Z)
        Ch = Z'*[-s.stage{n}.Xi(actX,1:2);
                 -C(:,1:2)];
        phih = Z'*[-s.stage{n}.xi(actX);-phi];
	else
        Ch = zeros(0,2);
        phih = [];
	end
    
    V = [eye(2);
         uK(:,:,n)];
    kh = [zeros(2,1);uk(:,:,n)];
    
    Pc = V'*(W+Ph)*V;
    qc = ((qh'+kh'*(W+Ph))*V)';
end
    fprintf('\n')