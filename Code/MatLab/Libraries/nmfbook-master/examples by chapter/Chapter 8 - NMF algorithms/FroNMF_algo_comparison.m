% This code can be used to generate Figures 8.3 and 8.4. 
% It compares NMF algorithms on three data sets: CBCL (dense), classic and
% tdt2 (sparse). 
clc; clear all; 

% Load data set, and set rank.  
% Choose beteween experience equal to 
% 1 --> CBCL, r=49
% 2 --> CBCL, r=10
% 3 --> Classic, r=30
% 4 --> TDT2, r=30
experience = 1;
if experience == 1
    load CBCL;
    r = 49;
elseif experience == 2
    load CBCL;
    r = 10;
elseif experience == 3
    load classic;
    r = 30;
elseif experience == 4
    load tdt2_top30
    r = 30;
end
[m,n] = size(X); 
numinit = 1;          % Number of initializations; in the book: 30 
options.timemax = 5;  % Time limit; in the book: 30 
options.beta0 = 0;    % Extrapolation between updates? No 
options.maxiter = Inf;% No maximum number of iterations 
options.display = 0;  % No display of error evolution
options.accuracy = 0; % Do not stop algorithms when the error stagnates
% Parameters to stop inner iterations 
options.alpha = 0.5;
options.delta = 0.1;
% Keep error at time steps in the vector tmf
nummarker = 40; % number of markers used for plotting the evolution of the error
tmf = 0:options.timemax/nummarker:options.timemax; 
ebest = +inf; % error of best solution computed 
for i = 1 : numinit
    fprintf('***** Initialization %2.0f of %2.0f *****\n', i, numinit); 
    % Random init  
    H = rand(r,n);
    W = rand(m,r);
    options.init.W = W;
    options.init.H = H;
    % Run the different algorithms
    disp('Running MU...');
    options.algo = 'MUUP';
    options.inneriter = 1;
    [Wmi,Hmi,emi,tmi] = FroNMF(X,r,options);
    options.inneriter = 100;
    ebest = min(ebest, min(emi));
    if i == 1 
        tm = tmf; 
        em = sumte(emi,tmi,tmf); 
        em = em/numinit; 
    else
        emf = sumte(emi,tmi,tmf); 
        em = em + emf/numinit; 
    end
    disp('Running A-MU...');
    options.algo = 'MUUP';
    [Wmai,Hmai,emai,tmai] = FroNMF(X,r,options);
    ebest = min(ebest, min(emai));
    if i == 1 
        tma = tmf; 
        [ema] = sumte(emai,tmai,tmf); 
        ema = ema/numinit; 
    else
        emaf = sumte(emai,tmai,tmf); 
        ema = ema + emaf/numinit; 
    end
    disp('Running ANLS...');
    options.algo = 'ASET';
    [Wai,Hai,eai,tai] = FroNMF(X,r,options);
    ebest = min(ebest, min(eai));
    if i == 1 
        ta = tmf; 
        [ea] = sumte(eai,tai,tmf); 
        ea = ea/numinit; 
    else
        eaf = sumte(eai,tai,tmf); 
        ea = ea + eaf/numinit; 
    end
    disp('Running ALS...');
    options.algo = 'ALSH';
    [Wlsi,Hlsi,elsi,tlsi] = FroNMF(X,r,options);
    ebest = min(ebest, min(elsi));
    if i == 1 
        tls = tmf; 
        [els] = sumte(elsi,tlsi,tmf); 
        els = els/numinit; 
    else
        elsf = sumte(elsi,tlsi,tmf); 
        els = els + elsf/numinit; 
    end
    disp('Running A-HALS...');
    options.algo = 'HALS';
    [Whi,Hhi,ehi,thi] = FroNMF(X,r,options);
    ebest = min(ebest, min(ehi));
    if i == 1 
        th = tmf; 
        [eh] = sumte(ehi,thi,tmf); 
        eh = eh/numinit; 
    else
        ehf = sumte(ehi,thi,tmf); 
        eh = eh + ehf/numinit; 
    end
    disp('Running E-A-HALS...');
    options.beta0 = 0.5;
    [Whei,Hhei,ehei,thei] = FroNMF(X,r,options);
    ebest = min(ebest, min(ehei));
    options.beta0 = 0;
    if i == 1 
        the = tmf; 
        [ehe] = sumte(ehei,thei,tmf); 
        ehe = ehe/numinit; 
    else
        ehef = sumte(ehei,thei,tmf); 
        ehe = ehe + ehef/numinit; 
    end
    disp('Running FPGM...');
    options.algo = 'FPGM';
    [Wfi,Hfi,efi,tfi] = FroNMF(X,r,options);
    ebest = min(ebest, min(efi));
    if i == 1 
        tf = tmf; 
        [ef] = sumte(efi,tfi,tmf); 
        ef = ef/numinit; 
    else
        eff = sumte(efi,tfi,tmf); 
        ef = ef + eff/numinit; 
    end
    disp('Running ADMM...');
    options.algo = 'ADMM';
    options.delta = 0.01;
    [Wadi,Hadi,eadi,tadi] = FroNMF(X,r,options);
    options.delta = 0.1;
    ebest = min(ebest, min(eadi));
    if i == 1 
        tad = tmf; 
        [ead] = sumte(eadi,tadi,tmf); 
        ead = ead/numinit; 
    else
        eadf = sumte(eadi,tadi,tmf); 
        ead = ead + eadf/numinit; 
    end
end
disp('Done.'); 
scriptfigNMFalgo