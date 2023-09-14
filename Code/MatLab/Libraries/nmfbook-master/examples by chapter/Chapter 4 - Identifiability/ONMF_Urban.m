% ONMF on the hyperspectral Urban image  
clear all; clc; 
load('./data sets/Urban.mat'); 
[m,n] = size(X); 
r = 6; 
% SPA initialization
options.display = 0;
K = SPA(X,r,options); 
% Run EM-ONMF
options.display = 1;
options.init = X(:,K); 
[W,H,e] = alternatingONMF(X,r,options); 
% Display results
affichage(H([2 1 6 3 5 4],:)',3,307,307); 
figure; 
plot( e ); 
xlabel('Iterations','Interpreter','latex'); 
ylabel('Relative error: $||X-WH||_F / ||X||_F$','Interpreter','latex'); 
title('Orthogonal NMF on Urban','Interpreter','latex'); 