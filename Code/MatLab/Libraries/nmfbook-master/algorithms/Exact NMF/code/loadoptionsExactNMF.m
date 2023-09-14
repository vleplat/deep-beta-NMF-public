% Default parameters

options.heuristic   = 'rbr';      % Heuristic algorithm: 'ms1','ms2','rbr','sa','hybrid'
options.maxiter     = 10;         % Number of attempts
options.tolerance   = 1e-6;       % Precision target for the objective function (relative error)
options.stop        = 1;          % Stop if optimal solution found
options.display     = 'on';       % Display: 'on', 'off'
options.algonmf     = 'HALS';     % NMF algorithm to use: see FroNMF.m: 'HALS','MUUP','ANLS','ADMM','ALSH','FPGM'
options.rndtype     = 'sparse10'; % Type of random initialization: 'rndcube','sparse00','sparse10','sparse01','sparse11'

% Parameters for LocalImprovement
options.alpha          = 0.99;    % If the error has not decreased by a factor alpha, the function stops
options.convergiter    = Inf;     % Number of iterations of an NMF algorithm between each test 
options.delta_t        = 1;       % Time of execution of an NMF algorithm between each test 
options.timemaxconverg = 600;     % The function stops after 'timemaxconverg' seconds

% MS2 options
options.Kms2 = 200;               % Number of different random initializations
options.Nms2 = 20;                % Number of iterations of an NMF algorithm for each initialization

% RBR options
options.Krbr = 10;                % Number of different random initializations
options.Nrbr = 50;                % Number of iterations of an NMF algorithm for each attempt

% SA options
options.T0      = -1;             % Initial temperature: logspace(T0,Tf,beta)
options.Tf      = -4;             % Final temperature
options.beta    = 20;             % Number of different levels of temperature
options.Ksa     = 50;             % Number of attempts at each level of temperature
options.Nsa     = 100;            % Number of iterations of an NMF algorithm for each attempt
options.J       = 2;              % Number of columns/rows of the current solution (W,H) to perturb