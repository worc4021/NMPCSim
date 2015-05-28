function sol = minWrapper(P,q,C,phi,actU,s)

K = s.R+P(3,3);
L = [s.stage{s.m}.Xi(actU,3);
     C(:,3)];

nkappa = numel(find(actU));
 
LHSconst = [-q(3);
            s.stage{s.m}.xi(actU);
            -phi];
LHSvar = [-P(3,1:2);
          -s.stage{s.m}.Xi(actU,1:2);
          -C(:,1:2)];

%       [Q,R] = qr(L');
%         if and(~isempty(L),size(L,1)+1<=size(Q,2))
%             EW = eig(Q(:,size(L,1)+1:end)'*K*Q(:,size(L,1)+1:end));
%         elseif size(L,1)+1>size(Q,2)
%             EW = zeros(size(K));
%         else
%             EW = eig(K);
%         end
%         
%         fprintf(' ')
%         for i = 1:length(EW)
%             fprintf('%6.1f ',EW(i))
%         end
%         fprintf(':')
        
[RHSconst,RHSvar,RHSorth] = KKTSolver(K,L,LHSconst,LHSvar);


sol.uK = RHSvar(1,:);
sol.uk = RHSconst(1);


if nkappa>0
    sol.kappaK = RHSvar(2:2+nkappa-1,:);
    sol.kappak = RHSconst(2:2+nkappa-1);
    if ~isempty(RHSorth)
        sol.kappaZ = RHSorth(2:2+nkappa-1,:);
    else
        sol.kappaZ = [];
    end
    if ~isempty(C)
        sol.rhoK = RHSvar(2+nkappa:end,:);
        sol.rhok = RHSconst(2+nkappa:end);
        if ~isempty(RHSorth)
            sol.rhoZ = RHSorth(2+nkappa:end,:); 
        else
            sol.rhoZ = [];
        end
    else
        sol.rhoK = zeros(0,2);
        sol.rhok = [];
        sol.rhoZ = []; 
    end
else
    sol.kappaK = zeros(0,2);
    sol.kappak = [];
    sol.kappaZ = [];
    
    if ~isempty(C)
        sol.rhoK = RHSvar(2:end,:);
        sol.rhok = RHSconst(2:end);
        if ~isempty(RHSorth)
            sol.rhoZ = RHSorth(2:end,:); 
        else
            sol.rhoZ = [];
        end
    else
        sol.rhoK = zeros(0,2);
        sol.rhok = [];
        sol.rhoZ = []; 
    end
    
end
        