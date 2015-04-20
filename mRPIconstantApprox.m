clear all
close all
clc

parameters;


[~,A,B] = a2DfunNew(x0,u0);
[~,K,~,gammaSQ,~] = terminalController(A,B,eye(2),diag([2,1]),1);
PSI = A+B*K;


F = [1;-1]/maxU;


[~,fval] = fmincon(@(x)nonlinearitiesWrapper(x,1),[x0;u0]+.9*[-maxX;maxU],[],[],[],[],...
    [x0;u0]-[maxX;maxU],[x0;u0]+[maxX;maxU],[],optimoptions('fmincon','GradObj','on'));

wMin = fval;

[~,fval] = fmincon(@(x)nonlinearitiesWrapper(x,-1),[x0;u0]+.9*[-maxX;maxU],[],[],[],[],...
    [x0;u0]-[maxX;maxU],[x0;u0]+[maxX;maxU],[],optimoptions('fmincon','GradObj','on'));


wMax = -fval;

G = [eye(2);-eye(2)];
g = [maxW;wMax;maxW;-wMin];

opt = sdpsettings('solver','cplex','cachesolvers',1,'verbose',0);

Lambda = [F*K,-ones(2,1);
          0,0,-1];

lambda = ones(3,1);


LambdaNext = [];
lambdaNext = [];

iter = 1;
iterMax = 50;
while and(~altIsContained(Lambda,lambda,LambdaNext,lambdaNext),iter<iterMax)
    
    
    if iter~=1
        Lambda = LambdaNext;
        lambda = lambdaNext;
    end
    
    LambdaNext=zeros(size(Lambda));
    lambdaNext=zeros(size(lambda));
    
        for j = 1:length(lambda)
            curIdx(j) = true;
            [wAst, fval] = cplexlp(-Lambda(curIdx,1:2), G, g);
            
            LambdaNext(j,:) = Lambda(j,:)*blkdiag(PSI,1);
            lambdaNext(j) = lambda(j)-Lambda(curIdx,1:2)*wAst(:);
            curIdx(j) = false;
        end
    [LambdaNext,lambdaNext] = rowReduce([Lambda;LambdaNext],[lambda;lambdaNext]);
    iter = iter+1;
end

pLambda = LambdaNext;
plambda = lambdaNext;

load('mRPIset.mat');

result = zeros(size(plambda));
MINI = false(size(plambda));

for i = 1:length(plambda)
    [x, fval] = cplexlp(-pLambda(i,1:2), Lambda, lambda);
    res = pLambda(i,1:2)*x;
    if and(pLambda(i,3)==0,res<=plambda(i))
        MINI(i)=true;
    elseif pLambda(i,3)~=0
        aStar=(plambda(i)-res)/pLambda(i,3);
        if pLambda(i,3)<0
            MINI(i)=true;
            result(i)=aStar;
        else
            MINI(i)=true;
            result(i)=-aStar;
        end
    end
end

alphaStar = max(result);

% X = sdpvar(3,1);
% C = YSet(X,[Lambda*X(1:2)<=lambda,X(3)>=alphaStar,X(3)<=alphaStar*1.05],opt);
[V,v] = rowReduce(pLambda(:,1:2),plambda-pLambda(:,3)*alphaStar);
% P = YSet(X,[V*X(1:2)<=v,X(3)==alphaStar],opt);
% Pt = Polyhedron([pLambda;0,0,1],[plambda;alphaStar*1.05]);
% plot(Pt,'alpha',.3,C,'alpha',.3)
plot(Polyhedron(V,v),'alpha',.3,Polyhedron(Lambda,lambda),'alpha',.3)

[~, minAlpha] = cplexlp([0,0,1], pLambda, plambda);