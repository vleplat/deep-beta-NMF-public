% Iteration of heuristic SA:
% For each level of temperature and for each of the K attemps,
% the solution (W,H) is perturbed.
% The perturbed solution is refined and is accepted 
% with some probability depending on the error. 
function [W,H,e]=iteration_sa(X,r,options,W,H)   
    [m,n]=size(X);
    
    %Initialisation
    if nargin <= 3
        [W,H]=initialization(m,n,options,r);
    end

    nX   = norm(X,'fro');
    e    = norm(X-W*H,'fro')/nX;
    emin = e;
    Wmin = W;
    Hmin = H;

    T0        = options.T0;
    Tf        = options.Tf;
    beta      = options.beta;
    K         = options.Ksa;
    N         = options.Nsa;
    J         = options.J;
    tol       = options.tolerance;
    timelimit = Inf;
    algonmf   = options.algonmf;
    
    temperatures = logspace(T0,Tf,beta);
    
    % For each level of temperature
    for j=1:beta
        T=temperatures(j);
        
        % At each temperature, there are K attempts
        for i=1:K
            Wt = W;
            Ht = H;
            
            % Perturbation of the current solution: (W,H)->(Wt,Ht)
            vecrandperm = randperm(r);
            indices     = vecrandperm(1:min(r,J));
            [Wt(:,indices),Ht(indices,:)]=initialization(m,n,options,length(indices));
            
            optionsNMFalgo.timemax = timelimit;
            optionsNMFalgo.maxiter = N;
            init.W = Wt; init.H = Ht;
            optionsNMFalgo.init = init;
            optionsNMFalgo.algo = algonmf;
            optionsNMFalgo.display = 0; 
            [Wt,Ht] = FroNMF(X,size(Wt,2),optionsNMFalgo);
            %[Wt,Ht] = algoNMF(X,Wt,Ht,N,timelimit,algonmf);
            
            % Estimation of the quality of the perturbation
            et    = norm(X-Wt*Ht,'fro')/nX;
            delta = et-e;
            
            % With probability or if the perturbed solution is better, we keep it.
            if delta < 0 || rand < exp(-delta/T)
                e = et;
                W = Wt;
                H = Ht;
                if e<emin
                    emin = e;
                    Wmin = W;
                    Hmin = H;
                end
                if emin<tol
                    break;
                end
            end
        end
        if emin<tol
            break;
        end
    end
    e = emin;
    W = Wmin;
    H = Hmin;
end