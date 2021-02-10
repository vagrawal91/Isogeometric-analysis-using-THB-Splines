%% RunMe file for 2D Plots and Refinement

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

%% Plot 2D basis

% for k = 1
%     for l = objU.levelBas{k}.activeIndex
%         plot2DTHB(objU,k,l)
%         hold on;
%     end
%     xlabel('x')
%     ylabel('y')
%     zlabel('value');
% end

%% Plot one single basis function

% figure
% plotOne2DTHB(objU,1,2)
% hold on
% objFun2 = bSplBasFun(2,objU.levelBas{1});
% %objFun2.plotOneBasisFun(objFun2.generOneBasisFun)
% C2 = objFun2.generOneBasisFun;
% scatter3(objFun2.plotVector,zeros(objFun2.sP,1),C2,'k.')
% scatter3(zeros(objFun2.sP,1),objFun2.plotVector,C2,'r.')
% xlabel('x')
% ylabel('y')
% zlabel('value');

%% Refine and plot 2D basis again

% refArea = [1/2 3/2]; 
% objU.ThbRefinement1DML(1,refArea);
% refArea = [1/4 3/4]; 
% objU.ThbRefinement1DML(2,refArea);

% for k = 1:objU.level
%     for l = objU.levelBas{k}.activeIndex
%         PlotLevel(objU,k,l)
%         hold on;
%     end
%     xlabel('x')
%     ylabel('y')
%     zlabel('value');
% end
