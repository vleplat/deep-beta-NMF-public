% sparse NMF on CBCL; see Figure 5.5 

clear all; close all; clc;

load('CBCL.mat');
[m,n] = size(X);
r = 49;
options.timemax = Inf;
options.maxiter = 100;
% SPA init
options.display = 0; 
K = SPA(X',r,options); 
H0 = X(K,:); 
W0 = NNLS(X',H0'); 
options.W = W0;
options.H = H0;
disp('***    NMF without sparsity constraints     ***')
options.display = 1; 
[W,H,e,t] = sparseNMF(X,r,options);
% Start from the NMF solution 
disp('*** Projected sparse NMF with sparsity 0.85 ***')
options.W = W;
options.H = H;
options.sW = 0.85;
[Ws,Hs,es,ts] = sparseNMF(X,r,options);
% Display the basis images
affichage(W,7,19,19); title('NMF');
affichage(Ws,7,19,19); title('sparse NMF (0.85)');