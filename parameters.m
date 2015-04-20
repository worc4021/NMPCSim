Ts = 3e-2;
C = 1;
g = 9.81;
m = .1;

x0 = [.05;0];
u0 = sqrt(g*m/C)*x0(1);

maxX = [10e-3;.105];
maxU = 1e-2;
maxW = 1e-5;

Q = diag([2,1]);
R = 1;