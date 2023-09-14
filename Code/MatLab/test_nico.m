%% --------------
%%% Test on TDT2
%%%--------------
close all; clear all; clc; 

load TDT2_500words 
r = [10; 5]; % ranks of the deep factorizations
rngsee = 2023; % control random seed
X = X(:,K); 
X=full(X);
% X = X/10000;
X = X/10000;
X = X'; 
[m,n] = size(X); 
smalltest = 1:10:n; 
X = X(:,smalltest); 

disp('Running multilayer KL-NMF'); 
maxiit = 200; % number of iterations 
options.maxiter = maxiit; 
options.rngseed = rngsee; % control the random seed 
[W,H,e] = multilayerKLNMF(X',r,options); 

disp('Running deep KL-NMF');
options.outerit = maxiit/2; % half iterations for deep KL-NMF
options.maxiter = maxiit/2; % half iterations for initialization with multilayer KL-NMF
options.min_vol = 1;        % activate minvol
options.delta = ones(length(r),1);          
options.epsi = 10^-8;      % can be reduced to 10^-3 - 10^-4 to speed up if needed
options.thres = 10^-5;      % stopping criterion for ADMM-procedure
options.innerloop = 100;    % inner loop for Step 1 of ADMM-procedure
options.maxIterADMM = 200; 
options.rho = 200; 
options.accADMM = 1; 
%options.lambda = [10; 1]; 
options.alpha_tilde = 1/2*ones(1,2); 
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


%%% sanity check for normalization constraints
if options.min_vol
    norm(sum(Wl{1},1) - ones(size(sum(Wl{1},1))))
    norm(sum(Wl{2},1) - ones(size(sum(Wl{2},1))))
else
    norm(sum(Hl{1},2) - ones(size(sum(Hl{1},2))))
    norm(sum(Hl{2},2) - ones(size(sum(Hl{2},2))))
end

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