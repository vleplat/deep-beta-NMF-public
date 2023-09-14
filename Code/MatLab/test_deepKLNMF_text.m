%% ------------------------------------------------------------------------
%%% Test on HSI data sets
%%%------------------------------------------------------------------------
close all; clear all; clc; 
cd('./');
addpath(genpath('./'));

%% ------------------------------------------------------------------------
%%% Libraries
%%%------------------------------------------------------------------------
% This code needs the toolboxes
% : from https://gitlab.com/ngillis/nmfbook/
% : from "DeepMF_flexibleFramework" framework
% -> both are copied within "Libraries" subdirectory 
% -> TIP : in main directory, type "addpath(genpath(pwd))" in Command Windo.


%% ------------------------------------------------------------------------
%%% Data set selection
%%%------------------------------------------------------------------------
%%% Urban data set
% load Urban
% X = A'; 
% mx = 307
% my = 307
% r = [6; 4];

% %%% small classic data set 
load('Datasets_TopicModel\classic400.mat')
%load('C:\Users\Nicolas\Dropbox\Data\text mining\ohscal.mat')
X = dtm'; 
subsamp = 500; 
K = SPA(X',subsamp); % subsample rows, could also be done in other ways
X = X(K,:); 
%[a,b] = sort(sum(X'),'descend'); 
%X = X(b(1:subsamp),:); 


% X = max(X,1e-6); 

r(2) = max(classid); 
r(1) = r(2)*2; 

%%% Tumor data set
% load('Tumor.mat') 
% X = M; 
% mx = 13;
% my = 11;
% r = [3; 2];

%%% Jasper Ridge data set
% load('jasperRidge2_R198.mat') 
% X = Y; 
% mx = 100;
% my = 100;
% r = [6; 4]; 

%% ------------------------------------------------------------------------
%%% Call of Algorithms
%%%------------------------------------------------------------------------
%%% main parameters
%rng(2023); 
options.maxiter = 200;  % max nu. of it. for init. stage
options.outerit = 100;  % max nu. of it. for our Algorithms
options.min_vol = 1;    % 0: Algorithm-1, 1: Algorithm-2
options.epsi = 10^-5;   % Algorithm 1: 10^-5, Algorithm 2: 10^-4 - 10^-3

%%% min-vol parameters
options.delta = ones(1,length(r));
options.alpha_tilde = 0.1*[1; 1];

% Tumor : options.alpha_tilde = [0.000005;0.01];
% Urban : options.alpha_tilde = [0.0000005;0.1];
% Moffet : options.alpha_tilde = [0.4;0.05]; [4;0.005]; [10;0.001];
% Jasper Ridge: options.alpha_tilde = [0.00004;0.01];

%%%% parameters for ADMM procedure
%%%%% Fast configuration - Algorithm 2 faster and accuracies of ADMM-procedure increase 
%%%%%  along iterations of the global scheme
options.rho = 100;           % 10-100: to fine tune
options.thres = 10^-4;      % stopping criterion for ADMM-procedure
options.innerloop = 1;      % inner loop for Step 1 of ADMM-procedure
options.maxIterADMM = 200;  
%%%%% Slow configuration - Algorithm 2 slower but accuracies of ADMM-procedure higher
%%%%% at each step 
% options.rho = 100;
% options.thres = 10^-4;
% options.innerloop = 3;
% options.maxIterADMM = 500;

%%% Call of methods
options.rngseed = 2000; 
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
    norm(sum(Hl{1},2) - ones(size(sum(Hl{1},2))))
    norm(sum(Hl{2},2) - ones(size(sum(Hl{2},2))))
end

%%% Errors
figure; 
plot(el./repmat(el(1,:)+1e-16,size(el,1),1))
title('Errors divided by errors at iteration 1'); 
legend('Level 1', 'Level 2', ...   %'Level 3', 
        'Weighted'); 

% classification accuracy of multilayer vs. deep 
H = inWH.H{2}*inWH.H{1}; 
[~,classDNMF] = max(H);
classDNMF = bestMap(classid,classDNMF); 
fprintf('Classfication accuracy of multilayer NMF is %2.0f%%\n', 100*sum(classDNMF == classid)/length(classid)); 

H = Hl{2}*Hl{1}; 
[~,classDNMF] = max(H);
classDNMF = bestMap(classid,classDNMF); 
fprintf('Classfication accuracy of deep NMF is %2.0f%%\n', 100*sum(classDNMF == classid)/length(classid)); 

% classification accuracy of NMF 
%[W,H,e] = betaNMF(X,r(2)); 
[Wl,Hl,el,inWH,output] = deepKL_NMF(X,r(2),options);
[~,classDNMF] = max(H);
classDNMF = bestMap(classid,classDNMF); 
fprintf('Classfication accuracy of NMF is %2.0f%%\n', 100*sum(classDNMF == classid)/length(classid)); 

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