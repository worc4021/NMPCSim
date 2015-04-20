function [dx, dfdx, dfdu] = a2DfunNew(x,u)

parameters;


dx = x + Ts*[x(2);g-C/m*(u/x(1))^2];
    
dfdx = eye(2) + Ts*[0,1;2*C/(m*x(1)^3)*u^2,0];
    
dfdu = Ts*[0;-2*C*u/(m*x(1)^2)];