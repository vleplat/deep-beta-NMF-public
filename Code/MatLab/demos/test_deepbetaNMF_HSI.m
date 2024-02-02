%% ------------------------------------------------------------------------
%%% Test on HSI data sets
%%%------------------------------------------------------------------------
close all; clear all; clc; 

%% ------------------------------------------------------------------------
%%% Libraries
%%%------------------------------------------------------------------------
% This code needs the toolboxes
% : from https://gitlab.com/ngillis/nmfbook/
% -> It is copied within "Libraries" subdirectory 
% -> TIP : in main directory, type "addpath(genpath(pwd))" in Command Windo.


%% ------------------------------------------------------------------------
%%% Data set selection
%%%------------------------------------------------------------------------
%%% Urban data set
% load Urban
% % X = A'; 
% mx = 307;
% my = 307;
% r = [6 2];
%%%X = X/max(X(:));

% %%% Moffet data set
load Moffet
mx = 50;
my = 50;
r = [4 2];

%%% Tumor data set
% load('Tumor.mat') 
% X = M; 
% mx = 13;
% my = 11;
% r = [3 2];

%%% Jasper Ridge data set
% load('jasperRidge2_R198.mat') 
% X = Y; 
% mx = 100;
% my = 100;
% r = [6 4]; 

%%% Samson data set
% load('samson_1.mat') 
% X = V; 
% mx = 95;
% my = 95;
% r = [4 2];
% X=max(X,10^-9);

%% ------------------------------------------------------------------------
%%% Call of Algorithms
%%%------------------------------------------------------------------------
%%% main parameters
rng(2023); 
options.beta = 3/2;
options.HnormType = 'cols'; % "cols": for e^T H = e^T, "rows": for H e = e 
options.maxiter = 1000;  % max nu. of it. for init. stage
options.outerit = 1000;  % max nu. of it. for our Algorithms
options.min_vol = 0;    % 0: Algorithm-1, 1: Algorithm-2
options.epsi = 10^-7;   % Algorithm 1: 10^-5, Algorithm 2: 10^-4 - 10^-3

%%% min-vol parameters
options.delta = ones(1,length(r));
options.alpha_tilde = [4;1];
% Tumor : [2.5*10^-5;10^-2];
% Urban : [0.0000005;0.1];
% Moffet : [4;1];
% Jasper Ridge: [0.00004;0.01];
% Samson: [0.45;0.45]; or [1.50;0.15];

%%%% parameters for ADMM procedure
options.rho = 100;           % 10-100: to fine tune
options.thres = 10^-6;      % stopping criterion for ADMM-procedure
options.innerloop = 1;      % inner loop for Step 1 of ADMM-procedure
options.maxIterADMM = 200;  
options.accADMM = 1;        % Call for the accelerated ADMM procedure

%%% Call of methods
[Wl,Hl,el,inWH,output] = deepKL_NMF(X,r,options);

%% ------------------------------------------------------------------------
%%% Post-processing
%%%------------------------------------------------------------------------
close all
%%% sanity check for normalization constraints
if options.min_vol
    norm(sum(Wl{1},1) - ones(size(sum(Wl{1},1))))
    norm(sum(Wl{2},1) - ones(size(sum(Wl{2},1))))
else
    if strcmp(options.HnormType,'rows')
        norm(sum(Hl{1},2) - ones(size(sum(Hl{1},2))))
        norm(sum(Hl{2},2) - ones(size(sum(Hl{2},2))))
    elseif strcmp(options.HnormType,'cols')
        norm(sum(Hl{1},1) - ones(size(sum(Hl{1},1))))
        norm(sum(Hl{2},1) - ones(size(sum(Hl{2},1))))
    end
end

%%% Errors
figure; 
plot(el./repmat(el(1,:)+1e-16,size(el,1),1))
title('Errors divided by errors at iteration 1'); 
legend('Level 1', 'Level 2', ...   %'Level 3', 
        'Weighted'); 

%%% Display abundance maps 
affichage(Hl{1}',2,mx,my); title('First layer deep KL-NMF'); 
affichage((Hl{2}*Hl{1})',2,mx,my); title('Second layer deep KL-NMF'); 
% affichage((Hl{3}*Hl{2}*Hl{1})',2,mx,my); title('Third layer deep KL-NMF'); 

affichage(inWH.H{1}',2,mx,my);  title('First layer multilayer KL-NMF');  
affichage((inWH.H{2}*inWH.H{1})',2,mx,my); title('Second layer multilayer KL-NMF');  
% affichage((inWH.H{3}*inWH.H{2}*inWH.H{1})',2,mx,my); title('Third layer multilayer KL-NMF');  

%%% Spectral signatures
figure;
plot(Wl{1});
title('First layer deep KL-NMF'); 
figure;
plot(Wl{2});
title('Second layer deep KL-NMF'); 

%%% min-vol info
if options.min_vol
    figure;
    semilogy(output.e_m); hold on;
    semilogy(output.logdetEvol);
    legend('$f(W,H)$','$\log \det(W_1^TW_1 + \delta I)$','$\log \det(W_2^TW_2 + \delta I)$','Interpreter','latex')
    grid on;

    fprintf(' ->Final betadivergence: %0.2f \n', el(end));
    fprintf(' ->Final penalty term: %0.2f \n', output.e_m(end));
    fprintf(' ->Initial ratio : %0.2f \n', output.ratio(1,:));
    fprintf(' ->Final ratio : %0.2f \n', output.ratio(2,:));

end