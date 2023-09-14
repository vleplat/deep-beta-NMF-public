% Compare min-vol models on the Urban hyperspectral image; 
% see Example 4.48 in the book 
clear all; clc; 

load Urban; 
r = 6; 

options.lambda = 1;
options.maxiter = 1000; % In the book: 1000
options.target = 0.05; 

% min-vol NMF (1) 
options.model=1; 
disp('Running min-vol NMF (1)...'); 
[W1,H1,e1,er11,er21] = minvolNMF(X,r,options); 

% min-vol NMF (2)  
options.model=2; 
disp('Running min-vol NMF (2)...'); 
[W2,H2,e2,er12,er22] = minvolNMF(X,r,options); 

% min-vol NMF (3) 
options.model=3; 
disp('Running min-vol NMF (3)...'); 
[W3,H3,e3,er13,er23] = minvolNMF(X,r,options); 


% Display extracted abundance maps, that is, rows of W
affichage(H1',6,307,307); title('min-vol NMF (1)'); 
affichage(H2',6,307,307); title('min-vol NMF (2)'); 
affichage(H3',6,307,307); title('min-vol NMF (3)'); 

% Display extracted spectral signatures, that is, columns of W
figure; 
subplot(1,3,1); 
plot(W1); 
title('min-vol NMF (1)', 'Interpreter','Latex'); 
leg1 = legend('1','2','3','4','5','6'); 
set(leg1,'Interpreter','latex'); 
subplot(1,3,2); 
plot(W2); 
title('min-vol NMF (2)', 'Interpreter','Latex'); 
leg2 = legend('1','2','3','4','5','6'); 
set(leg2,'Interpreter','latex'); 
subplot(1,3,3); 
plot(W3); 
title('min-vol NMF (3)', 'Interpreter','Latex'); 
leg3 = legend('1','2','3','4','5','6'); 
set(leg3,'Interpreter','latex'); 