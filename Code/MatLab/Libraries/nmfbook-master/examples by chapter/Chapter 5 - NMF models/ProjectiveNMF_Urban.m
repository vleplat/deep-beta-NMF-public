% Projective NMF on Urban; see Figure 5.2
clear all; clc; 
load Urban; 
[m,n] = size(X); 
r = 6; 
optionsSPA.display = 0; 
K = SPA(X',r,optionsSPA); 
H0 = X(K,:)'; 
options.init = H0; 
options.maxiter = 500; 
[H,e] = projectiveNMF(X',r,options); 
affichage(H,3,307,307); 
title('Solution of projective NMF'); 
figure; 
plot(e); 
xlabel('Iterations','Interpreter','latex'); 
ylabel('$\| X - WW^\top X \|_F / \|X\|_F$','Interpreter','latex'); 
title('projective NMF'); 