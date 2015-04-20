function sol = maxWrapper(P,q,Ch,phih,s)

% [K  L'][x]        =   [LHSconst]  + [LHSvar][theta]
% [L  0 ][lambda]   =   [LHSconst]  + [LHSvar][theta]

actW = s.wIdx(:,s.m);

[C,~] = chol(-P(1:2,1:2)+s.gammasq*eye(2));
fprintf('%9.2f ',diag(-C));
fprintf(': ');

nlambda = length(find(actW));

K = blkdiag(-s.gammasq*eye(2),P(1:2,1:2));
L = [-s.D,eye(2);
     s.G(actW,:),zeros(nlambda,2);
     zeros(size(Ch,1),2),Ch(:,1:2)];
 
LHSconst = [zeros(2,1);
            -q;
            zeros(2,1);
            s.g(actW);
            phih];
            
LHSvar = [zeros(2,3);
          zeros(2,3);
          s.A,s.B;
          s.H(actW,:);
          zeros(size(Ch,1),3)];
          
%         [Q,R] = qr(L');
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

sol.wM = RHSvar(1:2,:);
sol.wm = RHSconst(1:2);

sol.xpM = RHSvar(3:4,:);
sol.xpm = RHSconst(3:4);

sol.muM = RHSvar(5:6,:);
sol.mum = RHSconst(5:6);

if ~isempty(Ch)
    sol.muZ = RHSorth(5:6,:);
else
    sol.muZ = [];
end

if nlambda>0
    sol.lambdaM = RHSvar(7:7+nlambda-1,:);
    sol.lambdam = RHSconst(7:7+nlambda-1);

    
    if ~isempty(RHSorth)
        sol.lambdaZ = RHSorth(7:7+nlambda-1,:);
    else
        sol.lambdaZ = [];
    end
    
    if ~isempty(Ch)
    sol.zetaM = RHSvar(7+nlambda+1:end,:);
    sol.zetam = RHSconst(7+nlambda+1:end);
        if ~isempty(RHSorth)
            sol.zetaZ = RHSorth(7+nlambda+1:end,:);
        else
            sol.zetaZ = [];
        end
    else
        sol.zetaM = [];
        sol.zetam = [];
        sol.zetaZ = [];
    end
else
   sol.lambdaM = zeros(0,3);
   sol.lambdam = zeros(0,1);
   sol.lambdaZ = zeros(0,0);
   
   if ~isempty(Ch)
    sol.zetaM = RHSvar(7:end,:);
    sol.zetam = RHSconst(7:end);
        if ~isempty(RHSorth)
            sol.zetaZ = RHSorth(7:end,:);
        else
            sol.zetaZ = [];
        end
    else
        sol.zetaM = [];
        sol.zetam = [];
        sol.zetaZ = [];
    end
end