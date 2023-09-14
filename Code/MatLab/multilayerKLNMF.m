% Deep KL-NMF

%% * Input *
% X : m-by-n nonnegative data matrix
% r : vector of ranks
% Options: See betaNMF.m

%% * Output *
% W,H are arrays, that is,
% X \approx W{1}*H{1}, W{1} \approx W{2}*H{2}, ..., W{L-1} \approx W{L}*H{L} are the factorizations of the
% L levels.
% el(i) is the error at level i, and el(L+1) is the global error

function [W,H,e] = multilayerKLNMF(X,r,options)

if nargin <= 2
    options = [];
end
if ~isfield(options,'maxiter')
    options.maxiter = 500;
end
if ~isfield(options,'display')
    options.display = 1;
end
if ~isfield(options,'normalize')
    options.normalize = 2; %2: scale (W,H) such that He = e
                           %3: scale (W,H) such that W^Te = e
end
options.beta = 1;
L = length(r); 
options.delta = 1e-6; 

% Initialization by multiplayer KL-NMF
for i = 1 : L
    if i == 1
        %% SNPA-based initialization 
        %[Ksnpa,Hsnpa] = SNPA(X,r(i),options); 
        %options.W = full(X(:,Ksnpa)); 
        %options.H = Hsnpa;
        if isfield(options,'rngseed')
            rng(options.rngseed); 
        end
        [W{1},H{1},ei] = betaNMF(X,r(i),options);
        [W{1},H{1}] = normalizeWH(W{1},H{1},options.normalize); 
        e(1,i) = ei(end);
    else
        %[Ksnpa,Hsnpa] = SNPA(W{i-1},r(i),options); 
        %options.W = W{i-1}(:,Ksnpa); 
        %options.H = Hsnpa;
        if isfield(options,'rngseed')
            rng(options.rngseed); 
        end
        [W{i},H{i},ei] = betaNMF(W{i-1},r(i),options);
        [W{i},H{i}] = normalizeWH(W{i},H{i},options.normalize); 
        e(1,i) = ei(end);
    end
    if options.display == 1
        fprintf('Layer %2.0d done.\n', i);
    end
end