if ~exist('noRecompute','var')
    clear all
    close all
    clc
end

parameters;

[~,A,B] = a2DfunNew(x0,u0);


if 0
    n = 5;
    [PT,K,~,gammaSQ,~] = terminalController(A,B,eye(2),Q,R);
    PSI = A+B*K;
    
    [X1,U] = meshgrid(linspace(-maxX(1),maxX(1),n),linspace(-maxU,maxU,n));
    Hx1 = zeros(n^2,3);
    Hx2 = zeros(n^2,3);
    for i = 1:n
        for j = 1:n
                x = [X1(i,j);0];
                u = U(i,j);
                [~,Y1] = ode45(@(t,y)selectorWrapperNew(t,x,u,x0,u0,1),[0,1],zeros(4,1));
                [~,Y2] = ode45(@(t,y)selectorWrapperNew(t,x,u,x0,u0,2),[0,1],zeros(4,1));

                Hx = [Y1(end,1:2)',Y2(end,1:2)',Y1(end,3:4)']-[A,B];
                Hx2(n*(i-1)+j,:) = Hx(2,:);
        end
    end
    
    [P2,p2] = facetEnumeration(Hx2(:,[1,3]),ones(n^2,1));

    
    V2 = vertexReduction(Hx2,ones(n^2,1));
    
    save('data.mat','Hx2','P2','p2','V2','PSI','K','Q','R','PT','gammaSQ');
else
    load('data.mat');
end

G = [[1;-1],zeros(2,1);
    zeros(2,1),[1;-1];
    zeros(length(V2),1),ones(length(V2),1);
    zeros(length(V2),1),-ones(length(V2),1)];
g = [ones(4,1)*maxW;zeros(length(G)-4,1)];
H = [zeros(4,3);V2;-V2];

if ~exist('noRecompute','var')
    idx = (1:length(V2))';

    F = [1;-1]/maxU;

    Lambda = [[eye(2);-eye(2)]/diag(maxX);
              F*K];

    lambda = ones(6,1);

    LambdaNext = [];
    lambdaNext = [];

    iter = 1;
    iterMax = 50;
    while and(~altIsContained(Lambda,lambda,LambdaNext,lambdaNext),iter<iterMax)

        if iter~=1
            Lambda = LambdaNext;
            lambda = lambdaNext;
        end

        for i = 1:length(idx)
            curIdx = false(length(lambda),1);
            M = [zeros(1,3);V2(idx(i),:)]*[eye(2);K];
            for j = 1:length(lambda)
                curIdx(j) = true;
                [wAst, fval] = cplexlp(-Lambda(curIdx,:), G, g);
                tmp = (PSI+M);

                LambdaHelp = [Lambda(j,:)*PSI;      % both PWA are constant
                              Lambda(j,:)*[PSI(1,:);tmp(2,:)]];   % first constant, second linear
                lambdaHelp = [lambda(j)-Lambda(curIdx,:)*wAst(:);
                              lambda(j)-Lambda(curIdx,:)*[wAst(1);0]];
                [LambdaHelp,lambdaHelp] = inequalityReduction([Lambda;LambdaHelp],[lambda;lambdaHelp]);
                [LambdaNext,lambdaNext] = inequalityReduction([LambdaNext;LambdaHelp],[lambdaNext;lambdaHelp]);


                curIdx(j) = false;
            end
        end

        [LambdaNext,lambdaNext] = inequalityReduction([Lambda;LambdaNext],[lambda;lambdaNext]);
        iter = iter+1;
    end
%     figure(2)
%     plot(Polyhedron(LambdaNext,lambdaNext))

    Lambda = LambdaNext;
    lambda = lambdaNext;

    save('mRPIset.mat','Lambda','lambda')

end