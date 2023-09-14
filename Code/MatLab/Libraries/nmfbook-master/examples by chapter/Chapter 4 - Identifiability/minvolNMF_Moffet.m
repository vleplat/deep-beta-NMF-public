% Comparison of min-vol models on the Moffet hyperspectral image. 
%
% In this numerical experiment, we compare two min-vol NMF models: 
% One where H is normalized so that the entries in each column sum to at
% most one, that is, H^T e <= e, and the other where the entries in each 
% column sum one, that is, H^T e = e. 
% This is an interesting example because the Moffet data set contains
% pixels with very small spectral signatures; they correspond to water in
% the image. 
% The sum-to-at-most-one constraint considers this spectral signature as
% background noise, while using the sum-to-one constraint allows to
% identify the pixels containing water. 
% 
% For more information on this dataset; see 
% Dobigeon, N., Moussaoui, S., Coulon, M., Tourneret, J. Y., & Hero, A. O., 
% MCMC algorithms for supervised and unsupervised linear unmixing of 
% hyperspectral images. European Astronomical Society Publications Series, 
% 59, 381-401, 2013. 
% https://www.researchgate.net/publication/350549057_MCMC_Algorithms_for_Supervised_and_Unsupervised_Linear_Unmixing_of_Hyperspectral_Images

clear all; clc; 

load Moffet; 
r = 3; 

options.lambda = 1;
options.maxiter = 300; 
options.target = 0.05; 

% min-vol NMF (1) --> H^T e <= e 
options.model=1; 
disp('Running min-vol NMF with H^T e <= e ...'); 
[W1,H1,e1,er11,er21] = minvolNMF(X,r,options); 

% min-vol NMF (4) --> H^T e == e 
options.model=4; 
disp('Running min-vol NMF with H^T e  = e...'); 
[W4,H4,e4,er14,er24] = minvolNMF(X,r,options); 

% Display extracted abundance maps, that is, rows of H
affichage(H1',3,50,50); title('min-vol NMF - $H^T e \leq e$', 'Interpreter','Latex'); 
affichage(H4',3,50,50); title('min-vol NMF - $H^T e  = e$', 'Interpreter','Latex'); 

% Display extracted spectral signatures, that is, columns of W
figure; 
subplot(1,2,1); 
plot(W1); 
title('$H^T e \leq e$', 'Interpreter','Latex'); 
leg1 = legend('vegetation','soil','soil + vegetation'); 
set(leg1,'Interpreter','latex'); 
subplot(1,2,2); 
plot(W4); 
title('$H^T e  = e$', 'Interpreter','Latex'); 
leg2 = legend('water','soil','vegetation'); 
set(leg2,'Interpreter','latex'); 
