options.maxiter = 200;  % max nu. of it. for init. stage
options.outerit = 100;  % max nu. of it. for our Algorithms
options.min_vol = 1;    % 0: Algorithm-1, 1: Algorithm-2
options.epsi = 10^-5;   % Algorithm 1: 10^-5, Algorithm 2: 10^-4 - 10^-3
options.normalize = 2+options.min_vol; 

%%% min-vol parameters
options.delta = ones(1,length(r));
options.alpha_tilde = 0.001*[1; 1];

%%%% parameters for ADMM procedure
%%%%% Fast configuration - Algorithm 2 faster and accuracies of ADMM-procedure increase 
%%%%%  along iterations of the global scheme
options.rho = 100;           % 10-100: to fine tune
options.thres = 10^-4;      % stopping criterion for ADMM-procedure
options.innerloop = 1;      % inner loop for Step 1 of ADMM-procedure
options.maxIterADMM = 200;  