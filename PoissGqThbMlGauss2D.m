function PoissGqThbMlGauss2D(nEL,p)

clc;
a = -2;
b = 2;
N = nEL;
resol = 0.01;
lvl = 2;
objU = thbSplBasML(a,b,p,N,resol,lvl);
objV = thbSplBasML(a,b,p,N,resol,lvl);

f = @(x,y) 2; 
g = @(x,y) 10*x+(x^2+y^2)/2;

dBC = boundCond('Dirichlet','Dirichlet',0,0);

psU = PoissSolvNEW(objU,f);
psV = PoissSolvNEW(objV,f);
[Stiffn, rhs, ~, ~,~] = assembleMl(psU,psV);
sol = solveSyst(psU,psV,Stiffn,rhs,dBC);
uhU = generSolThb(psU,psV,sol);

% figure(1)
% [X,Y] = meshgrid(objU.levelBas{1}.knotVector,objV.levelBas{1}.knotVector);
% Z=ones(length(X),length(Y));
% s = mesh(X,Y,Z);

% figure(1)
% [X,Y] = meshgrid(objU.levelBas{1}.plotVector,objU.levelBas{1}.plotVector);
% Z=ones(length(X),length(Y)).*uhU;
% s = mesh(X,Y,Z);
% s.FaceColor = 'flat';
% plotBas(objU)

% refArea = [-1 1]; 
% objU.ThbRefinement1DML(1,refArea);
% objV.ThbRefinement1DML(1,refArea);
%              
% psU = PoissSolvNEW(objU,f);
% psV = PoissSolvNEW(objV,f);
% [Stiffn, rhs, ~, ~,~] = assembleMl(psU,psV);
% sol = solveSyst(psU,psV,Stiffn,rhs,dBC);
% uhU = generSolThb(psU,psV,sol);
% figure(2)
% [X,Y] = meshgrid(objU.levelBas{1}.plotVector,objU.levelBas{1}.plotVector);
% Z=ones(length(X),length(Y)).*uhU;
% s = mesh(X,Y,Z);
% s.FaceColor = 'flat';

% figure(2)
% [X,Y] = meshgrid(objU.levelBas{2}.knotVector,objV.levelBas{2}.knotVector);
% Z=ones(length(X),length(Y));
% s = mesh(X,Y,Z);
% plotBas(objU);

end