% Heuristic for Exact NMF. Given a nonnegative matrix X and a factorization
% rank r, it tries to find W>=0 with r columns and H>=0 with r rows such
% that X = WH. If it does not achieve this goal, it returns the solution
% found with the smallest error ||X-WH||_F. 
%
% See A. Vandaele, N. Gillis, F. Glineur and D. Tuyttens,
% "Heuristics for Exact Nonnegative Matrix Factorization", arXiv, 2014. 
% If you use the code, please cite the paper.
% See http://sites.google.com/site/nicolasgillis/ 
%
% [W,H,e,t] = exactNMFheur(X,r,options)
%
% Input.
%   X              : (m x n) matrix to factorize
%   r              : nonnegative rank
%   options        : optional field with parameters (see loadoptionsExactNMF.m file)
% Output.
%   (W,H)    : nonnegative matrices s.t. WH approximates M
%   (e,t)    : error and CPU time of each attempt (see options.maxiter)

function [W,H,e,t] = exactNMFheur(X,r,options)

    % The rank-one case can be solved via the SVD
    if r == 1
        disp('For r=1, we use the best rank-one approximation via the SVD.');
        f = cputime;
        [u,s,v] = svds(X,1);
        W = abs(u)*sqrt(s);
        H = abs(v')*sqrt(s);
        e = norm(X-W*H,'fro')/norm(X,'fro');
        t = cputime-f;
        fprintf('--> relative error = %2.2f%%\n', e*100);
        return;
    end

    if nargin<3
        options=update_options();
    else
        options=update_options(options);
    end
    if strcmp(options.display,'on')==1
        display(logdisplay_parameters(X,r,options));
    end
    
    emax=Inf;
    for i=1:options.maxiter
        f = cputime;
        switch options.heuristic
            case 'ms1'
                [Wi,Hi,ei,e0] = heuristic_ms1(X,r,options);
            case 'ms2'
                [Wi,Hi,ei,e0] = heuristic_ms2(X,r,options);
            case 'rbr'
                [Wi,Hi,ei,e0] = heuristic_rbr(X,r,options);
            case 'sa'
                [Wi,Hi,ei,e0] = heuristic_sa(X,r,options);
            case 'hybrid'
                [Wi,Hi,ei,e0] = heuristic_hybrid(X,r,options);
        end
        t(i) = cputime - f;
        e(i) = ei;
        if ei<emax
            emax  = ei;
            [W,H] = normNMF(Wi,Hi);
        end
        
        if strcmp(options.display,'on')==1
            display(logdisplay_attempts(i,ei,options));
        end
        
        if options.stop==1 && ei<options.tolerance
            break;
        end
    end
end

function options=update_options(optnew)
    if nargin==0
        loadoptionsExactNMF;
    else
        loadoptionsExactNMF;
        if isfield(optnew,'heuristic')
            options.heuristic=optnew.heuristic;
        end
        if isfield(optnew,'maxiter')
            options.maxiter=optnew.maxiter;
        end
        if isfield(optnew,'tolerance')
            options.tolerance=optnew.tolerance;
        end
        if isfield(optnew,'stop')
            options.stop=optnew.stop;
        end
        if isfield(optnew,'display')
            options.display=optnew.display;
        end
        if isfield(optnew,'algonmf')
            options.algonmf=optnew.algonmf;
        end
        if isfield(optnew,'rndtype')
            options.rndtype=optnew.rndtype;
        end
        
        if isfield(optnew,'alpha')
            options.alpha=optnew.alpha;
        end
        if isfield(optnew,'convergiter')
            options.convergiter=optnew.convergiter;
        end
        if isfield(optnew,'delta_t')
            options.delta_t=optnew.delta_t;
        end
        if isfield(optnew,'timemaxconverg')
            options.timemaxconverg=optnew.timemaxconverg;
        end
        
        if isfield(optnew,'Kms2')
            options.Kms2=optnew.Kms2;
        end
        if isfield(optnew,'Nms2')
            options.Nms2=optnew.Nms2;
        end
        
        if isfield(optnew,'Krbr')
            options.Krbr=optnew.Krbr;
        end
        if isfield(optnew,'Nrbr')
            options.Nrbr=optnew.Nrbr;
        end
        
        
        if isfield(optnew,'T0')
            options.T0=optnew.T0;
        end
        if isfield(optnew,'Tf')
            options.Tf=optnew.Tf;
        end
        if isfield(optnew,'beta')
            options.beta=optnew.beta;
        end
        if isfield(optnew,'Ksa')
            options.Ksa=optnew.Ksa;
        end
        if isfield(optnew,'Nsa')
            options.Nsa=optnew.Nsa;
        end
        if isfield(optnew,'J')
            options.J=optnew.J;
        end
    end
end

function s=logdisplay_parameters(X,r,options)
    s = sprintf('Factorizing a %dx%d matrix using r=%d with %d attempts of',size(X,1),size(X,2),r,options.maxiter);
    switch options.heuristic
        case 'ms1'
            s=sprintf('%s MS1\n',s);
            s=sprintf('%sParameters: ',s);
        case 'ms2'
            s=sprintf('%s MS2\n',s);
            s=sprintf('%sParameters: K=%d,N=%d - ',s,options.Kms2,options.Nms2);
        case 'rbr'
            s=sprintf('%s RBR\n',s);
            s=sprintf('%sParameters: K=%d,N=%d - ',s,options.Krbr,options.Nrbr);
        case 'sa'
            s=sprintf('%s SA\n',s);
            s=sprintf('%sParameters: K=%d,N=%d,T0=%d,Tf=%d,beta=%d,J=%d - ',s,options.Ksa,options.Nsa,options.T0,options.Tf,options.beta,options.J);
        case 'hybrid'
            s=sprintf('%s HYBRID\n',s);
            s=sprintf('%sParameters: Krbr=%d,Nrbr=%d,Ksa=%d,Nsa=%d,T0=%d,Tf=%d,beta=%d,J=%d - ',s,options.Krbr,options.Nrbr,options.Ksa,options.Nsa,options.T0,options.Tf,options.beta,options.J);
    end
    s=sprintf('%salgo=%d,rndtype=%s',s,options.algonmf,options.rndtype);
    s=sprintf('%s,LI(alpha=%1.2g,delta_t=%d)',s,options.alpha,options.delta_t);
end

function s=logdisplay_attempts(i,ei,options)    
    s = sprintf('\tAttempt %4d/%d',i,options.maxiter);
    if ei<options.tolerance
        s=sprintf('%s --> EXACT factorization with relative error = %1.3g%%',s,100*ei);
    else
        s=sprintf('%s --> relative error = %1.3g%%',s,100*ei);
    end
end
