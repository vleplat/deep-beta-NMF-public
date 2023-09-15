% Test on CBCL 
clear all; clc; 

addpath(genpath(pwd))
% This code also needs the toolbox from https://gitlab.com/ngillis/nmfbook/

load CBCL; 
% X = X/100; % for numerical stability 
r = [80; 40; 20] % ranks of the deep factorizations
rngsee = 35; % control random seed
disp('Running multilayer KL-NMF'); 
maxiit = 1000; % number of iterations
options.maxiter = maxiit; 
options.rngseed = rngsee; % control the random seed 
[W,H,e] = multilayerKLNMF(X',r,options); 

disp('Running deep KL-NMF');
options.outerit = maxiit/2; % half iterations for deep KL-NMF
options.maxiter = maxiit/2; % half iterations for initialization with multilayer KL-NMF
% for min-vol
options.min_vol = 0;        % activate minvol
options.delta = ones(length(r),1);          
options.epsi = 10^-10;      % can be reduced to 10^-3 - 10^-4 to speed up if needed
%options.alpha_tilde = 0.2; % ex:0.001 (none) - 0.5 (significant)
options.thres = 10^-4;      % stopping criterion for ADMM-procedure
options.innerloop = 100;    % inner loop for Step 1 of ADMM-procedure
options.maxIterADMM = 200;  
[Wl,Hl,el,inWH,output] = deepKL_NMF(X',r,options); 

%% Display results 
% Figure default
set(0, 'DefaultAxesFontSize', 25);
set(0, 'DefaultLineLineWidth', 2);
maxiy = 0; 
plot(el(:,1)/e(1)); hold on; 
maxiy = max(maxiy, max(el(:,1)/e(1))); 
plot(el(:,2)/e(2),'--'); 
maxiy = max(maxiy, max(el(:,2)/e(2))); 
plot(el(:,3)/e(3),'-.'); 
maxiy = max(maxiy, max(el(:,3)/e(3))); 
legend('Level 1', 'Level 2', 'Level 3','Interpreter','latex'); 
   % 'Weighted', 'ML-NMF-1', 'ML-NMF-2', 'ML-NMF-3'); 
xlabel('Iterations','Interpreter','latex'); 
ylabel('Ratio deep vs.\ multilayer','Interpreter','latex'); 
axis([0 options.outerit 0 ceil(10*max(el(:,1)/e(1)))/10])
grid on

numparligne = 10; 
affichage([Hl{1}' ones(size(X,1),numparligne+mod(r(1),10)) ...
           Hl{1}'*Hl{2}'  ones(size(X,1),numparligne+mod(r(2),10)) ... 
           Hl{1}'*Hl{2}'*Hl{3}'] , numparligne,19,19); 
title('Basis images extracted by deep KL-NMF'); 

affichage([H{1}' ones(size(X,1),numparligne+mod(r(1),10)) ...
           H{1}'*H{2}'  ones(size(X,1),numparligne+mod(r(2),10)) ... 
           H{1}'*H{2}'*H{3}'] , numparligne,19,19);
title('Basis images extracted by multilayer KL-NMF'); 

disp('Sparsity of multilayer KL-NMF:') 
[sp_col(H{1}'), sp_col(H{1}'*H{2}'), sp_col(H{1}'*H{2}'*H{3}')] 

disp('Sparsity of deep KL-NMF:')  
[sp_col(Hl{1}'), sp_col(Hl{1}'*Hl{2}'), sp_col(Hl{1}'*Hl{2}'*Hl{3}')] 
