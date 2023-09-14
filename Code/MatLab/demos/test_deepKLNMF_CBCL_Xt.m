% Test on CBCL 
clear all; clc; 

% addpath("Utils/") 
% addpath(genpath(pwd))
% This code also needs the toolbox from https://gitlab.com/ngillis/nmfbook/


load CBCL 
X = X/100;

rngsee = 1984; 

r = [80; 40; 20]
maxiit = 200; 
options.maxiter = maxiit; 
rng(rngsee); 
[W,H,e] = multilayerKLNMF(X',r,options)
options.outerit = maxiit/2*4; 
% for min-vol
options.min_vol = 0;         % activate minvol
options.delta = ones(length(r),1);          
options.epsi = 10^-10;        % can be reduced to 10^-3 - 10^-4 to speed up if needed
%options.alpha_tilde = 0.2;   % ex:0.001 (none) - 0.5 (significant)
options.thres = 10^-4;      % stopping criterion for ADMM-procedure
options.innerloop = 200;      % inner loop for Step 1 of ADMM-procedure
options.maxIterADMM = 200;  

rng(rngsee); 
[Wl,Hl,el,inWH,output] = deepKL_NMF(X',r,options); 

%% Post-processing
figure; 
plot(el./repmat(el(1,:),size(el,1),1))
title('Errors divided by errors at iteration 1'); 
hold on; 
plot([1 options.maxiter], [e(1)/el(1,1) e(1)/el(1,1)],'--'); 
plot([1 options.maxiter], [e(2)/el(1,2) e(2)/el(1,2)],'--'); 
plot([1 options.maxiter], [e(3)/el(1,3) e(3)/el(1,3)],'--'); 
legend('Level 1', 'Level 2', 'Level 3', ... 
    'Weighted', 'ML-NMF-1', 'ML-NMF-2', 'ML-NMF-3'); 

numparligne = 10; 
% affichage([Wl{1} ones(size(X,1),numparligne+mod(r(1),10)) ...
%            Wl{1}*Wl{2}  ones(size(X,1),numparligne+mod(r(2),10)) ... 
%            Wl{1}*Wl{2}*Wl{3}] , numparligne,19,19); 
affichage([Hl{1}' ones(size(X,1),numparligne+mod(r(1),10)) ...
           Hl{1}'*Hl{2}  ones(size(X,1),numparligne+mod(r(2),10)) ... 
           Hl{1}'*Hl{2}*Hl{3}] , numparligne,19,19); 
title('Basis images extracted by deep KL-NMF'); 

% affichage([W{1} ones(size(X,1),numparligne+mod(r(1),10)) ...
%            W{1}*W{2}  ones(size(X,1),numparligne+mod(r(2),10)) ... 
%            W{1}*W{2}*W{3}] , numparligne,19,19); 
affichage([H{1}' ones(size(X,1),numparligne+mod(r(1),10)) ...
           H{1}'*W{2}  ones(size(X,1),numparligne+mod(r(2),10)) ... 
           H{1}'*W{2}*W{3}] , numparligne,19,19);
title('Basis images extracted by multilayer KL-NMF'); 


% Sparsity of multilayer KL-NMF 
% [sp_col(W{1}), sp_col(W{1}*W{2}), sp_col(W{1}*W{2}*W{3})] 
[sp_col(H{1}'), sp_col(H{1}'*W{2}), sp_col(H{1}'*W{2}*W{3})] 


% Sparsity of multilayer KL-NMF 
% [sp_col(Wl{1}), sp_col(Wl{1}*Wl{2}), sp_col(Wl{1}*Wl{2}*Wl{3})] 
[sp_col(Hl{1}'), sp_col(Hl{1}'*Wl{2}), sp_col(Hl{1}'*Wl{2}*Wl{3})] 


% min-vol info
if options.min_vol
    figure;
    semilogy(output.e_m); hold on;
    semilogy(output.logdetEvol);
    legend('$f(W,H)$','$\log \det(W_1^TW_1 + \delta I)$','Interpreter','latex')
    grid on;

    fprintf(' ->Final betadivergence: %0.2f \n', el(end));
    fprintf(' ->Final penalty term: %0.2f \n', output.e_m(end));
    fprintf(' ->Initial ratio : %0.2f \n', output.ratio(1));
    fprintf(' ->Final ratio : %0.2f \n', output.ratio(2));

end