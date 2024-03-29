% Test on CBCL 
clear all; clc; 

addpath(genpath(pwd))
% This code also needs the toolbox from https://gitlab.com/ngillis/nmfbook/

load CBCL; 
options.beta = 3/2;
options.HnormType = 'cols'; % "cols": for e^T H = e^T, "rows": for H e = e 
options.epsi = 10^-7;
r = [80, 40, 20, 10] % ranks of the deep factorizations
rngsee = 35; % control random seed
disp('Running multilayer beta-NMF'); 
maxiit = 2000; % number of iterations
options.maxiter = maxiit; 
options.rngseed = rngsee; % control the random seed 
[W,H,e] = multilayerKLNMF(X',r,options); 

disp('Running deep beta-NMF');
options.outerit = maxiit/2; % half iterations for deep beta-NMF
options.maxiter = maxiit/2; % half iterations for initialization with multilayer beta-NMF
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
plot(el(:,4)/e(4),'-.'); 
maxiy = max(maxiy, max(el(:,4)/e(4))); 
legend('Level 1', 'Level 2', 'Level 3','Level 4','Interpreter','latex'); 
   % 'Weighted', 'ML-NMF-1', 'ML-NMF-2', 'ML-NMF-3'); 
xlabel('Iterations','Interpreter','latex'); 
ylabel('Ratio deep vs.\ multilayer','Interpreter','latex'); 
axis([0 options.outerit 0 ceil(10*max(el(:,1)/e(1)))/10])
grid on

numparligne = 10; 
affichage([Hl{1}' ones(size(X,1),numparligne+mod(r(1),10)) ...
           Hl{1}'*Hl{2}'  ones(size(X,1),numparligne+mod(r(2),10)) ... 
           Hl{1}'*Hl{2}'*Hl{3}' ones(size(X,1),numparligne+mod(r(3),10)) ... 
           Hl{1}'*Hl{2}'*Hl{3}'*Hl{4}'] , numparligne,19,19); 
title('Basis images extracted by deep beta-NMF'); 

affichage([H{1}' ones(size(X,1),numparligne+mod(r(1),10)) ...
           H{1}'*H{2}'  ones(size(X,1),numparligne+mod(r(2),10)) ... 
           H{1}'*H{2}'*H{3}' ones(size(X,1),numparligne+mod(r(3),10)) ...
           H{1}'*H{2}'*H{3}'*H{4}'] , numparligne,19,19);
title('Basis images extracted by multilayer beta-NMF'); 

disp('Sparsity of multilayer beta-NMF:') 
[sp_col(H{1}'), sp_col(H{1}'*H{2}'), sp_col(H{1}'*H{2}'*H{3}')] 

disp('Sparsity of deep beta-NMF:')  
[sp_col(Hl{1}'), sp_col(Hl{1}'*Hl{2}'), sp_col(Hl{1}'*Hl{2}'*Hl{3}')] 
