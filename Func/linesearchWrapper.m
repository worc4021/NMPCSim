function [s,XT,UT,WT] = linesearchWrapper(x0,xE,s)

tol = 1e-9;

d = xE-x0;
xC = x0;

XT = [];
UT = [];
WT = [];

t = 0;

iter = 0;
stepsize = 0;

fprintf('Iteration : Step size : State : Stage constraints : Disturbance constraints : Disturbance right hand side\n');
fprintf('%2.0d. : %9.2e : [ %8.4f , %8.4f ]  : ', iter,stepsize,xC(1),xC(2))
    for l = 1:3
        for m = 1:s.N
            if l == 1
                fprintf('%d',s.stage{m}.xiIdx')
                fprintf(' ')
            elseif l == 2
                fprintf('%d',s.wIdx(:,m)')
                fprintf(' ')
            elseif l==3
                fprintf('%d',s.rhsIdx(:,m)')
                fprintf(' ')
            end
        end
        fprintf(' : ')
    end
    iter = iter+1;


while d'*xC<d'*d-tol;

    if iter==18
        0;
    end
tic;
[wM,wm,xpM,xpm,lambdaM,lambdam,lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,...
    kappaK,kappak,kappaZ,rhoK,rhok,rhoZ] = solver(s);
TMP = toc;
if TMP>t
    t = TMP;
end

    [x,u,w,lambda,kappa] = trajectoryEvaluator(xC,wM,wm,xpM,xpm,lambdaM,...
        lambdam,lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,...
        kappaZ,rhoK,rhok,rhoZ,s);
    XT = [XT;x];
    UT = [UT;u];
    WT = [WT;w];
    
    
    [lSys,stepsize] = linesearch(xC,xE,wM,wm,xpM,...
    xpm,lambdaM,lambdam,lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,...
    kappaK,kappak,kappaZ,rhoK,rhok,rhoZ,s);


    if d'*(xC + stepsize*(xE-xC))>=xE'*xE-tol
        xC = xE;
        [x,u,w,lambda,kappa] = trajectoryEvaluator(xC,wM,wm,xpM,xpm,lambdaM,...
        lambdam,lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,...
        kappaZ,rhoK,rhok,rhoZ,s);
    XT = [XT;x];
    UT = [UT;u];
    WT = [WT;w];
        display('Terminated')
        break;
    end

    if and(stepsize<1,stepsize>0)
        s = lSys;
        s.t = t;
        xC = xC + stepsize*(xE-xC);
            fprintf('%2.0d. : %9.2e : [ %8.4f , %8.4f ]  : ', iter,stepsize,xC(1),xC(2))
            for l = 1:3
                for m = 1:s.N
                    if l == 1
                        fprintf('%d',s.stage{m}.xiIdx')
                        fprintf(' ')
                    elseif l == 2
                        fprintf('%d',s.wIdx(:,m)')
                        fprintf(' ')
                    elseif l==3
                        fprintf('%d',s.rhsIdx(:,m)')
                        fprintf(' ')
                    end
                end
                fprintf(' : ')
            end
            iter = iter+1;
    elseif stepsize>=1
        xC = xE;
        [x,u,w,lambda,kappa] = trajectoryEvaluator(xC,wM,wm,xpM,xpm,lambdaM,...
        lambdam,lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,...
        kappaZ,rhoK,rhok,rhoZ,s);
    XT = [XT;x];
    UT = [UT;u];
    WT = [WT;w];
        display('Terminated')
        break;
    else
        display('Infeasible problem')
        break;
    end

    
end
fprintf('\n')