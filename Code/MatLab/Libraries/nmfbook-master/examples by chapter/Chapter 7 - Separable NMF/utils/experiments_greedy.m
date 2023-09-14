% Comparison of greedy near-separable NMF algorithms
% 
%% *** Settings ***
% 1. Choose the type of synthetic data sets
% xp 
% 1 => well-conditioned `Dirichlet', 
% 2 => well-conditioned `Middle points', 
% 3 => ill-conditioned `Dirichlet', and 
% 4 => ill-conditioned `Middle points'. 
% condW = 4; % cond(W) = 10^condW

function [results, timings] = run_exp_greedy_sepNMF(m,n,r,xp,condW)

disp('*************************************************');
if xp == 1
    disp('    Experiment: well-conditioned Dirichlet');
elseif xp == 2
    disp('    Experiment: well-conditioned Middle points');
elseif xp == 3
    disp('    Experiment: ill-conditioned Dirichlet');
elseif xp == 4
    disp('    Experiment: ill-conditioned Middle points');
end
disp('***********************************************');
% 2. Noise levels 
numnoiselevels = 10; 
if xp <= 2
    delta=0:1/numnoiselevels:1; 
else
    %delta=0:1/numnoiselevels:1; 
    delta=logspace(-condW-2,0,numnoiselevels);  
end
% 3. Number of matrices generated per noise level 
nummat = 1; 
nalgo = 1:6; % tested algorithms
% Keep results 
results = zeros(length(nalgo), length(delta)); 
timings = zeros(length(nalgo), 1); 
% Do not display algorithms iterations
options.display = 0; 

%% *** Running the algorithms *** 
fprintf('Total number of noise levels: %2.0f \n', length(delta))
for i = 1 : length(delta)
    for k = 1 : nummat
        [W,H,Noise] = synthdatasetSepNMF(m,n,r,delta(i),xp,condW); 
        % Input matrix 
        Xt = W*H + Noise; 
        % Running the different algorithms
        for algoix = 1 : length(nalgo)
            e = cputime; 
            if nalgo(algoix) == 1 % SPA
                K = SPA(Xt,r,options); 
            elseif nalgo(algoix) == 2 % VCA
                [~, K,] = VCA(Xt,'Endmembers',r,'verbose','off'); 
            elseif nalgo(algoix) == 3 % FastAnchorWords 
                K = FastAnchorWords(Xt,r,options); 
            elseif nalgo(algoix) == 4 % SNPA
                options.maxitn = 500; 
                options.proj = 1; 
                K = SNPA(Xt,r,options); 
            elseif nalgo(algoix) == 5 % MVE-SPA
                K = MVESPA(Xt,r,options); 
            elseif nalgo(algoix) == 6 % SPA-SPA
                Kpre = SPA(Xt,r,options); 
                K = SPA(pinv(Xt(:,Kpre))*Xt,r,options); 
            end
            % Compute performance and running time 
            timings(algoix) = timings(algoix) + cputime-e; 
            rcur = sum(unique(K) <= r); 
            results(algoix,i) = results(algoix,i) + rcur/nummat/r; 
        end
    end
    fprintf('%1.0f...',i);
    if mod(i,10) == 0, fprintf('\n'); end
end