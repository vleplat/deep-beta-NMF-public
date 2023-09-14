% Projective NMF on CBCL; see Figure 5.3
clear all; clc; 
load CBCL;  
[m,n] = size(X); 
r = 49; 
K = SPA(X,r); 
W0 = X(:,K); 
options.W = W0; 
options.maxiter = 2000; 
[W,e] = projectiveNMF(X,r,options); 
affichage(W0,7,19,19); 
title('Initialization with separable NMF'); 
affichage(W,7,19,19); 
title('Solution of projective NMF'); 
figure; 
plot(e); 
xlabel('Iterations','Interpreter','latex'); 
ylabel('$\| X - WW^\top X \|_F$','Interpreter','latex'); 