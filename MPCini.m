clear all
close all
clc

noRecompute = 0;

the2DexpNew;
stageconstraints;

t = zeros(1,length(xiStage));
T = zeros(1,length(xiStage));

for N = 1:length(xiStage)

system.A = A;
system.B = B;
system.D = eye(2);
system.G = G;
system.g = g;
system.H = H;
system.Q = Q;
system.R = R;
system.gammasq = gammaSQ;
system.rhsIdx = false(length(g),N);
system.rhsIdx(1:3,:) = true(3,N);
system.rhsIdx(14,:) = true(1,N); % find(g)
system.wIdx = false(length(g),N);
system.P = PT;
system.N = N;
system.nuXi = 0;
system.xiS = 1;
system.stage = cell(1,N);

for m = (N-1):-1:0
    system.stage{m+1}.Xi = XiStage{end-m};
    system.stage{m+1}.xi = xiStage{end-m};
    system.stage{m+1}.xiIdx = false(size(xiStage{end-m}));
    system.nuXi = system.nuXi + length(xiStage{end-m});
    system.xiS = [system.xiS,system.nuXi];
end

xE = [2;50];
tic;
[outSys,XT,UT,WT] = linesearchWrapper([0;0],xE,system);
T(N) = toc;
t(N) = outSys.t;
end


figure(1)
hold on
for i = 0:(N-1)
    [VT,typeT] = vertexEnumeration(LambdaStage{end-i},lambdaStage{end-i});
    order = convhull(VT(:,1),VT(:,2));
    fill(VT(order,1),VT(order,2),'r')
    alpha(.2)
end

for i = 1:(size(XT,1)/2)
    plot(XT(2*i-1,:),XT(2*i,:),'b')
end
% 
% x0 = XT(end-1:end,1);
% xE = [0;-30];
% 
% [outSys1,XT1,UT1,WT1] = linesearchWrapper(x0,-xE,outSys);
% 
% for i = 1:(size(XT1,1)/2)
%     plot(XT1(2*i-1,:),XT1(2*i,:),'k')
% end
hold off
xlabel('$x_1$','interpreter','latex','fontsize',16);
ylabel('$x_2$','interpreter','latex','fontsize',16);


figure(2)
hold on
for i = 1:(N+1)
    plot(Polyhedron(blkdiag([1;-1],LambdaStage{end-i+1}),[(i-1)*[1;-1];lambdaStage{end-i+1}]),'alpha',.2)
end

for i = 1:(size(XT,1)/2)
    plot3(N:-1:0,XT(2*i-1,:),XT(2*i,:),'b')
end
hold off
xlabel('$m$','interpreter','latex','fontsize',16);
ylabel('$x_1$','interpreter','latex','fontsize',16);
zlabel('$x_2$','interpreter','latex','fontsize',16);
