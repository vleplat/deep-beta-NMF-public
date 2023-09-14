% Iteration of heuristic MS2:
% Generation of K different random initializations (W,H),
% N iterations of an NMF algorithm are applied to each pair.
% The output corresponds to the best pair.
function [W,H,e] = iteration_ms2(X,r,options)
    [m,n] = size(X);
    
    K         = options.Kms2;
    N         = options.Nms2;
    timelimit = Inf;
    algonmf   = options.algonmf;
    
    e  = Inf;
    W  = [];
    H  = [];
    nX = norm(X,'fro');
    
    for i=1:K
        [Wi,Hi] = initialization(m,n,options,r);
        
        optionsNMFalgo.timemax = timelimit; 
        optionsNMFalgo.maxiter = N; 
        init.W = Wi; init.H = Hi; 
        optionsNMFalgo.init = init; 
        optionsNMFalgo.algo = algonmf; 
        optionsNMFalgo.display = 0; 
        [Wi,Hi] = FroNMF(X,size(Wi,2),optionsNMFalgo); 
        %[Wi,Hi] = algoNMF(X,Wi,Hi,N,timelimit,algonmf);

        ei      = norm(X-Wi*Hi,'fro')/nX;
        if ei<e
            e = ei;
            W = Wi;
            H = Hi;
        end
    end
end