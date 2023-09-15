%% --------------
%%% Test on TDT2
%%%--------------
close all; clear all; clc; 
%% Keep most important words using NMF: this code comes from the NMF 
%% toolbox available from https://gitlab.com/ngillis/nmfbook/
load tdt2_top30 
X = X'; 
disp('*** Dataset tdt2_top30  ***'); 
fprintf('Sparsity: %2.2f%% of zero entries.\n', (1 - sum(X(:) > 0)/size(X,1)/size(X,2))*100) 
rng(2020); 
r = 20; 
[W,H] = FroNMF(X,r);
K = []; 
for i = 1 : r
    [a,b] = sort(-W(:,i));
    for j = 1 : 30 % Keep the 10 words with the largest value in W(:,i) 
        K = [K, b(j)]; 
    end
end 
K = unique(K); 
%% Parameters 
r = [20; 10; 5]; % ranks of the deep factorizations
rngsee = 2023; % control random seed
X = X(K,:); 
X=full(X);
X = X';  
disp('Running multilayer KL-NMF'); 
maxiit = 1000; % number of iterations 
options.maxiter = maxiit; 
options.rngseed = rngsee; % control the random seed 
[W,H,e] = multilayerKLNMF(X',r,options); 
disp('Running deep KL-NMF');
options.outerit = maxiit/2; % half iterations for deep KL-NMF
options.maxiter = maxiit/2; % half iterations for initialization with multilayer KL-NMF
options.min_vol = 0;        % activate minvol
options.epsi = 10^-9;       % can be reduced to 10^-3 - 10^-4 to speed up if needed
options.thres = 10^-4;      % stopping criterion for ADMM-procedure
options.innerloop = 100;    % inner loop for Step 1 of ADMM-procedure
options.maxIterADMM = 200; 
options.accADMM = 1; 
options.lambda = [4; 2; 1]; % weights the different layers 
[Wl,Hl,el,inWH,output] = deepKL_NMF(X',r,options); 
%% Display evolution of the error of deep KL-NMf vs. multilayer  
% Figure default
figure; 
set(0, 'DefaultAxesFontSize', 25);
set(0, 'DefaultLineLineWidth', 2);
maxiy = 0; 
for ri = 1 : length(r) 
    plotri = el(:,ri)/e(ri); 
    plot(plotri); hold on; 
    maxiy = max(maxiy, max(plotri));
end
legend('Level 1', 'Level 2', 'Level 3', 'Interpreter','latex'); 
xlabel('Iterations','Interpreter','latex'); 
ylabel('Ratio deep vs.\ multilayer','Interpreter','latex'); 
axis([0 options.outerit 0 ceil(10*maxiy)/10]); 
grid on

%% Construct Topics 
wordsK = words(K); 
for ri = 1 : length(r)
    Wc = W{ri};
    for i = 1 : r(ri)
        [a,b] = sort(Wc(:,i),'descend');
        for j = 1 : 10*ri % Keep the 10 words with the largest value in W(:,i)
            TopicsMulti{ri}{j,i} = wordsK{b(j)};
        end
    end
end
for ri = 1 : length(r)
    Wc = Wl{ri};
    for i = 1 : r(ri)
        [a,b] = sort(Wc(:,i),'descend');
        for j = 1 : 10*ri % Keep the 10 words with the largest value in W(:,i)
            TopicsDeep{ri}{j,i} = wordsK{b(j)};
        end
    end
end
disp('**************************************************************'); 
disp('You can vizualize the topics within TopicsMulti and TopicsDeep'); 
disp('**************************************************************'); 
