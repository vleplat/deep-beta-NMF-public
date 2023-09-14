% ONMF on CBCL 
clear all; clc; 
load('./data sets/CBCL.mat'); 
X = X';
[m,n] = size(X); 
r = 49; 
% SPA initialization
options.display = 0; 
K = SPA(X,49,options); 
% Run EM-ONMF
options.display = 1;
options.init = X(:,K); 
[W,H,e] = alternatingONMF(X,r,options);
% Display results: 
affichage(H',7,19,19); 
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 2);
figure; 
plot( e ); 
xlabel('Iterations'); 
ylabel('||X-WH||_F / ||X||_F'); 
title('Orthogonal NMF on CBCL','Interpreter','latex');  