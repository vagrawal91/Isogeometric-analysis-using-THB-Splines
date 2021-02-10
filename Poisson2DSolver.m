%% Poisson 2D Solver
% So far works only for N=1,2 case

clc;
clear all;

a = 0;
b = 2;
N = 1;
p = 2;
resol = 0.01;
lvl = 3;
objU = thbSplBasML(a,b,p,N,resol,lvl);
objV = thbSplBasML(a,b,p,N,resol,lvl);

% f = @(x,y) x^2+y^2-x-y; 
% g = @(x,y) x*(x-1)*y*(y-1);

f = @(x,y) 2; 
g = @(x,y) (1-x^2-y^2)/2;

dBC = boundCond('Dirichlet','Dirichlet',0,0);

psU = PoissSolv2D(objU,f);
psV = PoissSolv2D(objV,f);
[Stiffn, rhs, ~, ~,~] = assembleMl(psU,psV);
sol = solveSyst(psU,psV,Stiffn,rhs,dBC);
uhU = generSolThb(psU,psV,sol);
uhU2D = uhU*uhU';
%% Plot exact solution

% figure(1)
% [X,Y] = meshgrid(objU.levelBas{1}.knotVector,objV.levelBas{1}.knotVector);
% Z=ones(length(X),length(Y));
% s = surf(X,Y,(1-X^2-Y^2)/2);

%% Plot approximation

figure(2)
[X,Y] = meshgrid(objU.levelBas{1}.plotVector,objV.levelBas{1}.plotVector);
surf(X,Y,uhU2D);
% s.FaceColor = 'flat';
% pcolor(X,Y,Z)
%% Refine and plot again

refArea = [0 2]; 
objU.ThbRefinement1DML(1,refArea);
objV.ThbRefinement1DML(1,refArea);
psU = PoissSolv2D(objU,f);
psV = PoissSolv2D(objV,f);
[Stiffn, rhs, ~, ~,~] = assembleMl(psU,psV);
sol = solveSyst(psU,psV,Stiffn,rhs,dBC);
uhU = generSolThb(psU,psV,sol);
uhU2D = uhU*uhU';
figure(3)
[X,Y] = meshgrid(objU.levelBas{1}.plotVector,objV.levelBas{1}.plotVector);
surf(X,Y,uhU2D);
% plotBas(objU);

