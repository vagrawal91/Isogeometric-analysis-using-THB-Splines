%% Poisson 2D Solver
% Solver works perfectly, global refinement aswell, but local refinement 
% cannot be appplied in any case. For N=4,p=2 and refArea=[1 3/2] works
% fine though.

clc;
clear all;

a = 0;
b = 2;
N = 4;
p = 2;
resol = 0.01;
lvl = 3;
objU = thbSplBasML(a,b,p,N,resol,lvl);
objV = thbSplBasML(a,b,p,N,resol,lvl);

f = @(x,y) 2; 
g = @(x,y) (1-x^2-y^2)/2;
% f = @(x,y) 100*sin(pi*x); 

dBC = boundCond('Dirichlet','Dirichlet',0,0);

psU = PoissSolv2D(objU,f);
psV = PoissSolv2D(objV,f);
[Stiffn, rhs, ~, ~,~] = assembleMl(psU,psV);
sol = solveSyst(psU,psV,Stiffn,rhs,dBC);
uhU = generSolThb(psU,psV,sol);
uhU2D = uhU*uhU';

%% Plot approximation

% figure(1)
% [X,Y] = meshgrid(objU.levelBas{1}.plotVector,objV.levelBas{1}.plotVector);
% surf(X,Y,uhU2D);
% plotBas(objU);
% s.FaceColor = 'flat';
% figure(3)
% pcolor(X,Y,uhU2D)

%% Refine and plot again

refArea = [1 3/2]; 
objU.ThbRefinement1DML(1,refArea);
objV.ThbRefinement1DML(1,refArea);
psU = PoissSolv2D(objU,f);
psV = PoissSolv2D(objV,f);
[Stiffn, rhs, ~, ~,~] = assembleMl(psU,psV);
sol = solveSyst(psU,psV,Stiffn,rhs,dBC);
uhU = generSolThb(psU,psV,sol);
uhU2D = uhU*uhU';
figure(2)
[X,Y] = meshgrid(objU.levelBas{1}.plotVector,objV.levelBas{1}.plotVector);
surf(X,Y,uhU2D);
% figure(4)
% plotBas(objU);
% pcolor(X,Y,uhU2D)
