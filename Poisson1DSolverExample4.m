%% Poisson 1D Solver 3rd Example

clc;
clear all;

a = 0;
b = 2;
N = 6;
p = 2;
resol = 0.01;
lvl = 4;
obj = thbSplBasML(a,b,p,N,resol,lvl);

f = @(x) exp(x);
g = @(x) -exp(x)+7.3890561*x+1;
dBC = boundCond('Dirichlet','Neumann',0,0);

% plotBas(obj);
refArea = [1 2]; 
obj.ThbRefinement1DML(1,refArea);
% plotBas(obj);
refArea2 = [1.5 2];
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