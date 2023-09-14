% Iteration of heuristic RBR:
% There are r steps and at each step, we add a column/row to W/H.
% At the first step, the optimal solution can computed using the SVD.
function [W,H,e] = iteration_rbr(X,r,options)
    [m,n] = size(X);
    
    W=zeros(m,r);
    H=zeros(r,n);
    
    % For the rank-1 problem, the optimal solution is known.
    [w1,sig1,v1] = svds(X,1);
    W(:,1) = abs(w1);
    H(1,:) = sig1*abs(v1');

    % At each step, we add a column/row to W/H
    for i=2:r
        [W,H,e] = getRankPlusOne(i,X,W,H,options);
    end
end

% Adding a column and a row to respectively W and H:
% There are K attemps and the best one is kept.
function [W,H,e] = getRankPlusOne(ind,X,W,H,options)
    [m,n] = size(X);
    emin  = Inf;
    nX    = norm(X,'fro');
    
    K         = options.Krbr;
    N         = options.Nrbr;
    timelimit = Inf;
    algonmf   = options.algonmf;
    
    % There are K attempts
    for i=1:K
        Wt = W(:,1:ind);
        Ht = H(1:ind,:);
        
        % For each attempt, we randomly choose the last column and the last
        % row of respectively W and H.
        [Wt(:,ind),Ht(ind,:)]=initialization(m,n,options,1);

        optionsNMFalgo.timemax = timelimit; 
        optionsNMFalgo.maxiter = N; 
        init.W = Wt; init.H = Ht; 
        optionsNMFalgo.init = init; 
        optionsNMFalgo.algo = algonmf; 
        optionsNMFalgo.display = 0; 
        [Wt,Ht] = FroNMF(X,size(Wt,2),optionsNMFalgo); 
        
        e       = norm(X-Wt*Ht,'fro')/nX;
        
        if e<emin
            emin = e;
            Wmin = Wt;
            Hmin = Ht;
        end
    end
    e          = emin;
    W(:,1:ind) = Wmin;
    H(1:ind,:) = Hmin;
end