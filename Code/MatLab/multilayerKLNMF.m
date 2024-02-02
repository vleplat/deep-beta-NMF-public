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
if ~isfield(options,'beta')
    options.beta = 1;
end

if ~isfield(options,'HnormType')
    options.HnormType = 'rows';
end

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
        if strcmp(options.HnormType,'rows')
            [W{1},H{1},ei] = betaNMF(X,r(i),options);
            [W{1},H{1}] = normalizeWH(W{1},H{1},options.normalize); 
            e(1,i) = ei(end);
        elseif strcmp(options.HnormType,'cols')
            [m,n] = size(X); 
            W{1} = rand(m,r(i));
            H{1} = rand(r(i),n);
            % Display evolution of the iterations
            for k=1:options.maxiter
                prod = W{1}*H{1};
                Wt = W{1}';
                C = (Wt*(((prod).^(options.beta-2)).*X));
                D = Wt*(prod).^(options.beta-1);
                %%% Computation of Lagrange multipliers
                [~,I] = min(D,[],1);
                idx=sub2ind(size(D),I,1:n);
                H_current = H{1};
                mu_0_H = (D(idx)-C(idx).*(H_current(idx)))';
                mu_H = updatemu_hcols(C,D,H{1},1,mu_0_H,options.epsi);
                %%% Update matrix H_L 
                H{1} = H{1} .* (C./(D-ones(r(i),1)*mu_H'+eps));
                H{1} = max(H{1},eps);
                W{1} = MUbeta(X',H{1}',W{1}',options.beta)'; 

                if options.display == 1
                    if mod(i,10) == 0
                        fprintf('%1.2d...',i);
                    end
                    if mod(i,100) == 0
                        fprintf('\n');
                    end
                end
            end
            e(1,i) = betadiv(X,W{1}*H{1},options.beta);
        end
        
    else
        %[Ksnpa,Hsnpa] = SNPA(W{i-1},r(i),options); 
        %options.W = W{i-1}(:,Ksnpa); 
        %options.H = Hsnpa;
        if isfield(options,'rngseed')
            rng(options.rngseed); 
        end
        if strcmp(options.HnormType,'rows')
            [W{i},H{i},ei] = betaNMF(W{i-1},r(i),options);
            [W{i},H{i}] = normalizeWH(W{i},H{i},options.normalize); 
            e(1,i) = ei(end);
        elseif strcmp(options.HnormType,'cols')
            [m,n] = size(W{i-1}); 
            W{i} = rand(m,r(i));
            H{i} = rand(r(i),n);
            for k=1:options.maxiter
                prod = W{i}*H{i};
                Wt = W{i}';
                C = (Wt*(((prod).^(options.beta-2)).*W{i-1}));
                D = Wt*(prod).^(options.beta-1);
                %%% Computation of Lagrange multipliers
                [~,n] = size(W{i-1}); 
                [~,I] = min(D,[],1);
                idx=sub2ind(size(D),I,1:n);
                H_current = H{i};
                mu_0_H = (D(idx)-C(idx).*(H_current(idx)))';
                mu_H = updatemu_hcols(C,D,H{i},1,mu_0_H,options.epsi);
                %%% Update matrix H_L 
                H{i} = H{i} .* (C./(D-ones(r(i),1)*mu_H'+eps));
                H{i} = max(H{i},eps);
                W{i} = MUbeta(W{i-1}',H{i}',W{i}',options.beta)'; 
                if options.display == 1
                    if mod(i,10) == 0
                        fprintf('%1.2d...',i);
                    end
                    if mod(i,100) == 0
                        fprintf('\n');
                    end
                end
            end
            e(1,i) = betadiv(W{i-1},W{i}*H{i},options.beta);
        end
        
    end
    if options.display == 1
        fprintf('Layer %2.0d done.\n', i);
    end
end