% Deep KL-NMF

%% * Input *
% X : m-by-n nonnegative data matrix
% r : vector of ranks
% Options: See betaNMF.m

%% * Output *
% W and H are arrays, such that 
% X \approx W{1}*H{1}, W{1} \approx W{2}*H{2}, ..., W{L-1} \approx W{L}*H{L} 
%   are the factorizations of the L levels.
% el(i) is the error at level i, and el(L+1) is the global error

function [W,H,e,inWH,output] = deepKL_NMF(X,r,options)

if nargin <= 2
    options = [];
end
L = length(r); 
% Check that the ranks decrease 
[a,b] = sort(r,'descend'); 
if min(b == 1:L) == 0
    warning('The ranks of deep NMF should be decreasing.');
end
if ~isfield(options,'maxiter')
    options.maxiter = 500; % number of iterations for the init. 
end
if ~isfield(options,'outerit')
    options.outerit = 100;
end
if ~isfield(options,'display')
    options.display = 1;
end
if ~isfield(options,'min_vol')
    options.min_vol = 0;
    options.alpha = zeros(1,L);
end
if ~isfield(options,'delta')
    options.delta = ones(1,L);
end
if ~isfield(options,'epsi')
    options.epsi = 10^-6;
end
if ~isfield(options,'alpha_tilde')
    options.alpha_tilde = 0.05*ones(1,L);
end
if ~isfield(options,'accADMM')
    options.accADMM = 0;
end

% Misc parameters 
options.beta = 1;

if ~isfield(options,'W') || ~isfield(options,'H')
    % Initialization by multiplayer KL-NMF
    disp('Initialization with multi-layer KL-NMF');
    options.normalize = 2+options.min_vol; 
    [W,H,e] = multilayerKLNMF(X,r,options); 
end
%% Keep initial solutions in memory 
for i = 1 : L
    inWH.W{i} = W{i}; 
    inWH.H{i} = H{i}; 
end
%%  Misc. initializations
% Initialize the values of lambda, if not provided 
if ~isfield(options,'lambda')
    options.lambda = 1./e(1,:)'; 
end
% initial global error 
e(1,L+1) = options.lambda'*e(1,1:L)'; 

% Init. of quantities related to minvol approach
cumsum = options.lambda'*e(1,1:L)';
if options.min_vol 
    logdetSave = zeros(options.outerit+1,L); 
    output.ratio = zeros(2,L);
    for i=1:L
        % init. for alpha
        if i==1
            df = betadiv(X,W{i}*H{i},options.beta);
        else
            df = betadiv(W{i-1},W{i}*H{i},options.beta);
        end
        mv = log10(det(W{i}'*W{i}+options.delta(i)*eye(r(i))));
        options.alpha(i)=options.alpha_tilde(i)*df/abs(mv);

        % logdet and ratio of objective terms
        logdetSave(1,i) = mv;
        output.ratio(1,i) = options.alpha(i)*logdetSave(1,i)/df;

        % initial global error
        cumsum = cumsum+ options.alpha(i)*mv;  
    end 
    e_m(1,1) = cumsum;
    % initialization for mu (Lagrangian multipliers) - Layer L
    mu = zeros(r(L),1);
    % init. of matrix full for ones for updates of W_L
    E = ones(size(W{L}*H{L}));
end

%% Optimizing the levels alternatively
disp('Alternating optimization of the L layers in deep KL-NMF (iter : error)'); 
for it = 1 : options.outerit
    for i = 1 : L
        if i == 1 && L > 1
            if ~options.min_vol
                lam = options.lambda(2)/options.lambda(1); 
                [W{1},H{1}] = levelUpdateDeepKLNMF(H{1},X,W{1},W{2}*H{2},lam,options.epsi);
            else 
                Wp = W{2}*H{2};
                [W{1},H{1},res{i}] = levelUpdateDeepminvolKLNMF(H{1},X,W{1},Wp,options,i);
            end
            e(it+1,i) = betadiv(X,W{1}*H{1},options.beta);
            logdetSave(it+1,i) = log10(det(W{i}'*W{i}+options.delta(i)*eye(r(i))));
        elseif i == L
            if ~options.min_vol
                %%% Update of factors W_L and H_L
                H{i} = MUbeta(W{i-1},W{i},H{i},options.beta);
                W{i} = MUbeta(W{i-1}',H{i}',W{i}',options.beta)'; 
                %%% scale 
                [W{i},H{i}] = normalizeWH(W{i},H{i},2); 
                %%% error computation for the layer at hand
                e(it+1,i) = betadiv(W{i-1},W{i}*H{i},options.beta);
            else
                if i == 1
                    Wprev = X; 
                else
                    Wprev = W{i-1}; 
                end
                [m,n] = size(Wprev); 
                Y = inv(W{i}'*W{i} + options.delta(i)*eye(r(i)));
                Y_plus=max(0,Y);
                Y_minus=max(0,-Y);
                Phi=repmat(ones(1,n)*H{i}',m,1)-4*options.alpha(i)*W{i}*Y_minus;
                Eta=(8*options.alpha(i)*(W{i}*(Y_plus+Y_minus))).*((Wprev./(W{i}*H{i}+eps))*H{i}');
                Psi=4*options.alpha(i)*W{i}*(Y_plus+Y_minus);
                mu=updatemu(Phi,Eta,Psi,W{i},mu,options.epsi);
                %%% update maxtrix W ("dictionaries")
                W{i} = W{i} .*(((Phi+ones(m,1)*mu').^2+Eta).^(1/2)-(Phi+ones(m,1)*mu'))./(Psi+eps);
                W{i} = max(W{i},eps);
                %%% update maxtrix H 
                prod = W{i}*H{i};
                W_Lt = W{i}';
                H{i} = (H{i} .* (W_Lt*((E./prod).*Wprev))./(W_Lt*E+eps));
                H{i} = max(H{i},eps);
                %%% error computation for layer L 
                e(it+1,i) = betadiv(Wprev,W{i}*H{i},options.beta);
                logdetSave(it+1,i) = log10(det(W{i}'*W{i}+options.delta(i)*eye(r(i))));
            end
            
        else
            if ~options.min_vol
                lam = options.lambda(i+1)/options.lambda(i);
                [W{i},H{i}] = levelUpdateDeepKLNMF(H{i},W{i-1},W{i},W{i+1}*H{i+1},lam,options.epsi);
            else
                Wp = W{i+1}*H{i+1};
                [W{i},H{i},res{i}] = levelUpdateDeepminvolKLNMF(H{i},W{i-1},W{i},Wp,options,i);
            end
            e(it+1,i) = betadiv(W{i-1},W{i}*H{i},options.beta);
            logdetSave(it+1,i) = log10(det(W{i}'*W{i}+options.delta(i)*eye(r(i))));
        end % end of if for number of layer test
    end % end of loop on the layers
    e(it+1,L+1) = options.lambda'*e(it+1,1:L)';

    %errors relative to min-vol formulation
    if options.min_vol 
        e_m(it+1,1) = options.lambda'*e(it+1,1:L)'+ options.alpha*logdetSave(it+1,1:L)';%  global error
    end

    % Display evolution of the iterations
    if options.display == 1
        if mod(it,10) == 0
            fprintf('%2.0d : %1.2f - ',it,e(it+1,L+1));
            close all;
            if options.min_vol
                figure; semilogy(res{1});
            end
        end
        if mod(it,50) == 0
            fprintf('\n');
        end
    end
end % end of iterative loop

% Save outputs
output.e = e;

%errors relative to min-vol formulation
if options.min_vol 
    output.logdetEvol = logdetSave;
    output.e_m = e_m;
    output.alpha = options.alpha;
    for i=1:L
        df = e(end,i);
        output.ratio(2,i) = options.alpha(i)*logdetSave(end,i)/df;
    end
end
