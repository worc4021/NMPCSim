function [x,u,w,lambda,kappa] = trajectoryEvaluator(x0,wM,wm,xpM,xpm,lambdaM,lambdam,lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,kappaZ,rhoK,rhok,rhoZ,s)

x = zeros(2,s.N+1);
u = zeros(1,s.N);
w = zeros(2,s.N);
lambda = zeros(length(s.g),s.N);
kappa = cell(1,s.N);

x(:,1) = x0;

for i = 1:s.N
    u(:,i) = uK(:,:,i)*x(:,i)+uk(:,:,i);
    kappa{i} = kappaK{i}*x(:,i)+kappak{i};
    w(:,i) = wM(:,:,i)*[x(:,i);u(:,i)]+wm(:,:,i);
    lambda(:,i) = lambdaM(:,:,i)*[x(:,i);u(:,i)]+lambdam(:,:,i);
    x(:,i+1) = xpM(:,:,i)*[x(:,i);u(:,i)]+xpm(:,:,i);
    wMax = s.G*w(:,i)-diag(double(s.rhsIdx(:,i)))*(s.g+s.H*[x(:,i);u(i)]);
    xiMax = s.stage{i}.Xi*[x(:,i);u(i)]- s.stage{i}.xi;
end