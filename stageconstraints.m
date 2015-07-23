if ~exist('noRecompute','var')
    if ~exist('N','var')
        N = 50;
    end
    
    LambdaStage = cell(1,N+1);
    lambdaStage = cell(1,N+1);

    LambdaStage{end}=Lambda;
    lambdaStage{end}=lambda;

    XiStage = cell(1,N);
    xiStage = cell(1,N);

    Polyhedra = cell(1,N+1);
    Polyhedra{end} = Polyhedron(LambdaStage{end},lambdaStage{end});

    Xi = zeros(0,3);
    xi = zeros(0,1);

    SYS = [A,B];
    for k = N:-1:1
        for i = 1:length(idx)
            curIdx = false(length(lambdaStage{k+1}),1);
            M = [zeros(1,3);V2(idx(i),:)];
            for j = 1:length(lambdaStage{k+1})
                    curIdx(j) = true;
                    [wAst, fval] = cplexlp(-LambdaStage{k+1}(curIdx,:), G, g);
                    tmp = SYS + M;

                    XiHelp = [LambdaStage{k+1}(j,:)*SYS;      % both PWA are constant
                                  LambdaStage{k+1}(j,:)*[SYS(1,:);tmp(2,:)]];   % first constant, second linear
                    xiHelp = [lambdaStage{k+1}(j)-LambdaStage{k+1}(curIdx,:)*wAst(:);
                                  lambdaStage{k+1}(j)-LambdaStage{k+1}(curIdx,:)*[wAst(1);0]];
                    [Xi,xi] = inequalityReduction([Xi;XiHelp],[xi;xiHelp]);


                    curIdx(j) = false;
            end
        end
        [XiStage{k},xiStage{k}] = inequalityReduction([Xi;[zeros(2),F]],[xi;ones(2,1)]);
        Xi = zeros(0,3);
        xi = zeros(0,1);
        [LambdaStage{k},lambdaStage{k}] = projectPolyhedron(XiStage{k},xiStage{k},1);
        Polyhedra{k} = Polyhedron(LambdaStage{k},lambdaStage{k});
    end

    save('stageConstraints.mat','XiStage','xiStage','LambdaStage','lambdaStage');
    
    figure(1)
    hold on
    for i = 1:(N+1)
        plot(Polyhedra{i},'alpha',.3)
    end
    hold off
else
    load('stageConstraints.mat');
end