%% --------------
%%% Test on TDT2
%%%--------------
close all; clear all; clc; 

load TDT2_500words 
r = [10; 5]; % ranks of the deep factorizations
rngsee = 2023; % control random seed
X = X(:,K); 
X=full(X);
X = X'; 
[m,n] = size(X); 
smalltest = 1:10:n; 
X = X(:,smalltest); 

disp('Running multilayer KL-NMF'); 
maxiit = 300; % number of iterations 
options.maxiter = maxiit; 
options.rngseed = rngsee; % control the random seed 
[W,H,e] = multilayerKLNMF(X',r,options); 

disp('Running deep KL-NMF');
options.outerit = maxiit/2; % half iterations for deep KL-NMF
options.maxiter = maxiit/2; % half iterations for initialization with multilayer KL-NMF
options.min_vol = 0;        % activate minvol
options.delta = ones(length(r),1);          
options.epsi = 10^-9;      % can be reduced to 10^-3 - 10^-4 to speed up if needed
options.thres = 10^-4;      % stopping criterion for ADMM-procedure
options.innerloop = 100;    % inner loop for Step 1 of ADMM-procedure
options.maxIterADMM = 200; 
options.rho = 200; 
options.accADMM = 1; 
%options.lambda = [10; 1]; 
options.alpha_tilde = 0.01*ones(1,2); 
[Wl,Hl,el,inWH,output] = deepKL_NMF(X',r,options); 

%% Display results 
% Figure default
figure; 
set(0, 'DefaultAxesFontSize', 25);
set(0, 'DefaultLineLineWidth', 2);
firtiterdis = 1; 
maxiy = 0; 
plot(el(firtiterdis:end,1)/e(1)); hold on; 
maxiy = max(maxiy, max(el(firtiterdis:end,1)/e(1))); 
plot(el(firtiterdis:end,2)/e(2),'--'); 
maxiy = max(maxiy, max(el(firtiterdis:end,2)/e(2))); 
legend('Level 1', 'Level 2', 'Interpreter','latex'); 
xlabel('Iterations','Interpreter','latex'); 
ylabel('Ratio deep vs.\ multilayer','Interpreter','latex'); 
axis([0 options.outerit 0 ceil(10*maxiy)/10]); 
grid on

% %% ------------------------------------------------------------------------
% %%% Post-processing
% %%%------------------------------------------------------------------------
% % classification accuracy of multilayer vs. deep 
% [~,classDNMF] = max(W{1}');
% classDNMF = bestMap(classid,classDNMF); 
% fprintf('Classfication accuracy of multilayer NMF is %2.0f%%\n', 100*sum(classDNMF == classid)/length(classid)); 
% H = Hl{2}*Hl{1}; 
% [~,classDNMF] = max(H);
% classDNMF = bestMap(classid,classDNMF); 
% fprintf('Classfication accuracy of deep NMF is %2.0f%%\n', 100*sum(classDNMF == classid)/length(classid)); 
% % classification accuracy of NMF 
% [~,classDNMF] = max(Hl{1});
% classDNMF = bestMap(classid,classDNMF); 
% fprintf('Classfication accuracy of NMF is %2.0f%%\n', 100*sum(classDNMF == classid)/length(classid)); 