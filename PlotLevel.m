function PlotLevel(obj,lvl,ind)
        [X,Y] = meshgrid(obj.levelBas{lvl}.knotVector(1):obj.levelBas{lvl}.resol...
            :obj.levelBas{lvl}.knotVector(end),obj.levelBas{lvl}.knotVector(1)...
            :obj.levelBas{lvl}.resol:obj.levelBas{lvl}.knotVector(end));
        objFun1 = thbSplBasFun(ind,obj,lvl);
        CU = objFun1.generOneBasisFun();
        objFun1. plotOneBasisFun(CU)
        for k = obj.levelBas{lvl}.activeIndex
            objFun2 = bSplBasFun(k,obj.levelBas{lvl});
            CV = objFun2.generOneBasisFun();
           Z = CU(:,1)*CV(:,1)'; % tensor product structure
           surf(Y,X,Z);
           surf(X,Y,Z);
           hold on;
        end
end