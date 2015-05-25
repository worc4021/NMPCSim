clear all
close all
clc

p = path();
if isempty(strfind(p,'Func'));
    path(p,'Func');
end

noRecompute = 0;

the2DexpNew;
stageconstraints;

N = length(xiStage);

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
system.rhsIdx(1:4,:) = true(4,N);
system.wIdx = false(length(g),N);
system.P = PT;
system.N = N;
system.nuXi = 0;
system.xiS = 1;
system.stage = cell(1,N);

for m = N:-1:1
    system.stage{m}.Xi = XiStage{m};
    system.stage{m}.xi = xiStage{m};
    system.stage{m}.xiIdx = false(size(xiStage{m}));
    system.nuXi = system.nuXi + length(xiStage{m});
    system.xiS = [system.xiS,system.nuXi];
end

xE = [1;-50];

tic;
[outSys,XT,UT,WT] = linesearchWrapper([0;0],xE,system);
toc
figure(1)
hold on
for i = 1:N
    [VT,typeT] = vertexEnumeration(LambdaStage{i},lambdaStage{i});
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