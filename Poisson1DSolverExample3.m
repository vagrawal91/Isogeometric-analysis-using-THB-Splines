%% Poisson 1D Solver 3rd Example

clc;
clear all;

a = -1;
b = 1;
N = 8;
p = 1;
resol = 0.01;
lvl = 6;
obj = thbSplBasML(a,b,p,N,resol,lvl);

f = @(x) -exp(x);
g = @(x) exp(x)-1.1752*x-1.5431;
dBC = boundCond('Dirichlet','Dirichlet',0,0);

% plotBas(obj);
refArea = [-1/2 1]; 
obj.ThbRefinement1DML(1,refArea);
% plotBas(obj);
refArea2 = [-1/4 1/2];
obj.ThbRefinement1DML(2,refArea2);
% plotBas(obj);

ps = PoissSolv(obj,f);
[Stiffn, rhs, ~,~] = ps.assembleMl();
y = ps.solveSyst(Stiffn,rhs,dBC);
fplot(g,[obj.levelBas{1}.a obj.levelBas{1}.b],'b')
hold on;
uh =  ps.generSolThb(y);
plot(obj.levelBas{1}.plotVector(1:end-1),uh(1:end-1),'r')
hold on;
plot(obj.getAllKnots,0.01,'k*', 'markers',4)

% Errors in H1, L2 and Energy norms
[H1err,L2err,Enerr] = ps.errCalc(g,y);