function inte = selectorWrapperNew(t,x,u,x0,u0,selector)

Du = zeros(2);

[~,Dx,Du(:,1)] = a2DfunNew(x0+t*x,u0+t*u);

inte = [Dx(:,selector);Du(:,selector)];