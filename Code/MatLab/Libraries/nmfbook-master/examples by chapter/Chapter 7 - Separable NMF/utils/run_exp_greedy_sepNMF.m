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
% nalgo: tested algorithms 

function [results, timings] = run_exp_greedy_sepNMF(m,n,r,xp,condW,delta,nummat,nalgo,diri)

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
% Keep results 
results = zeros(length(nalgo), length(delta)); 
timings = zeros(length(nalgo), 1); 
% Do not display algorithms iterations
options.display = 0; 
options.maxitn = 500; % number of iterations within SNPA
options.dimred = 0; % MVE-SPA uses SPA as a preprocessing
%% *** Running the algorithms *** 
fprintf('Total number of noise levels: %2.0f \n', length(delta))
for i = 1 : length(delta)
    for k = 1 : nummat
        [W,H,Noise] = synthdatasetSepNMF(m,n,r,delta(i),xp,condW,diri); 
        % Input matrix 
        Xt = W*H + Noise; 
        % Permutation 
        permu = randperm(size(Xt,2)); 
        Xt = Xt(:,permu); 
        Kstar = find(permu <= r); 
        % Running the different algorithms
        ldr = 1;
        for algoix = 1 : length(nalgo)
            e = cputime; 
            if nalgo(algoix) == 1 % SPA
                K = SPA(Xt,r,options); 
            elseif nalgo(algoix) == 2 % VCA
                [~, K,] = VCA(Xt,'Endmembers',r,'verbose','off'); 
            elseif nalgo(algoix) == 3 % FastAnchorWords 
                K = FastAnchorWords(Xt,r,options); 
            elseif nalgo(algoix) == 4 % SNPA
                options.proj = 1; 
                K = SNPA(Xt,r,options); 
                Ksnpa = K; 
            elseif nalgo(algoix) == 5 % MVE-SPA
                if ldr == 1
                    options.dimred = 0;
                elseif ldr == 2
                    options.dimred = 1;
                elseif ldr == 3
                    options.dimred = 3;   
                end
                K = MVESPA(Xt,r,options); 
                ldr = ldr + 1; 
            elseif nalgo(algoix) == 6 % SPA-SPA
                Kpre = SPA(Xt,r,options); 
                K = SPA(pinv(Xt(:,Kpre))*Xt,r,options); 
            elseif nalgo(algoix) == 7 % FGNSR 
                [Y, Kf] = fgnsr(Xt, r, ones(size(Xt,2),1), 'maxiter', 1000, ... 
                    'verbose', false, 'proj', 'parproj', ... 
                    'col_weights', ones(size(Xt,2),1), 'dynamic_mu', 10); 
                K = SPA(Y',r,options); 
            end
            % Compute performance and running time 
            timings(algoix) = timings(algoix) + cputime-e; 
            rcur = length( intersect(K,Kstar) ); 
            results(algoix,i) = results(algoix,i) + rcur/nummat/r; 
        end
    end
    fprintf('%1.0f...',i);
    if mod(i,10) == 0, fprintf('\n'); end
end