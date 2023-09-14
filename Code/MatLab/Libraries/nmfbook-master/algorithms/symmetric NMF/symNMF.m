% Exact coordinate method for symmetric NMF (symNMF): given an n-by-n 
% symmetric nonnegative matrix A and a factorization rank r, it computes a
% n-by-r nonnegative matrix H, solving the problem
% 
%     min_{H >= 0}  1/2 * ||A - HH^T||_F^2
% 
% using an exact coordinate descent method. 
%
% See A. Vandaele, N. Gillis, Q. Lei, K. Zhong and I.S. Dhillon, "Efficient 
% and Non-Convex Coordinate Descent for Symmetric Nonnegative Matrix 
% Factorization", IEEE Trans. on Signal Processing 64 (21), pp. 5571-5584, 
% 2016. 
% 
% If you use the code, please cite the paper.
% The code is avaialble from https://sites.google.com/site/nicolasgillis/ 
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% To use symNMF, you need to compile the file symNMFmex.cpp using the 
% command: 
%              "mex -largeArrayDims symNMFmex.cpp"
% 
% For this, you need to install a compiler. 
% See, e.g., http://nl.mathworks.com/support/compilers/R2012a/win64.html 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
% [H,e,t] = symNMF(A,Hr,options)
%
% Input.
%   A              : (n x n) symmetric matrix to factorize
%   Hr             : inner rank of factorization OR initial matrix of size (n x r)
%   options        : optional field with parameters (see loadoptionsSymNMF.m file)
% Output.
%   H              : nonnegative matrices s.t. HH^T approximates A
%   (e,t)          : error 1/2 * ||A - HH^T||_F^2, and CPU time

function [H,e,t]=symNMF(A,Hr,options)
    e=[]; t=[];
    
    % A must be a symmetric matrix
    n = size(A,1);
    
    % Loading options
    if nargin<3
        options=update_options();
    else
        options=update_options(options);
    end
    
    run_possible=true;
    
    % Input argument #2 can be 'r' or the initial matrix 'H'
    if size(Hr,1)==1 && size(Hr,2)==1
        r = Hr;
        switch options.initmatrix
            case 'zeros'
                H = zeros(n,r);
            case 'dense01'
                if options.seed ~= -1
                    rng(options.seed);
                end
                H = rand(size(A,1),r);
            otherwise
                display(sprintf('Error - options.initmatrix is not correct (zeros or dense01).'));
                run_possible = false;
        end
    else
        H = Hr;
        r = size(H,2);
        if size(H,1)~=n
            display('Error - The size of A and the size of H must me the same');
            run_possible=false;
        end
    end
    
    %%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%
    % Maximum number of iterations
    if ~isnumeric(options.maxiter) || options.maxiter<=0 || floor(options.maxiter)-options.maxiter>0
        run_possible = false;
    else
        if options.maxiter>1e6
            display(sprintf('\nWarning: maxiter is > 1e6, we change it to 1e6.'));
            maxiter = 1e6;
        else
            maxiter = options.maxiter;
        end
    end
    
    % Maximum time of execution
    if ~isnumeric(options.timelimit) || options.timelimit<=0
        run_possible = false;
    else
        timelimit = options.timelimit;
    end
    
    % Shuffling the columns
    switch options.shuffle_columns
        case 0
            sc = 0;
        case 1
            sc = 1;
        otherwise
            display(sprintf('Error - options.sort_columns is not correct (0 or 1).'));
            run_possible = false;
    end
    if(run_possible)
        %%%%%%%%%%%% Scaling %%%%%%%%%%%%%%
        if sum(sum(H))>0
            nHtH    = norm(H'*H,'fro')^2;
            HtAHt   = sum(sum((H'*A).*(H')));
            scaling = HtAHt/nHtH;
            H       = sqrt(scaling)*H;
        end

        %%%%%%%%%%%% Initial objective function %%%%%%%%%%%%%%
        nA    = norm(A,'fro')^2;
        nHtH  = norm(H'*H,'fro')^2;
        HtAHt = sum(sum((H'*A).*(H')));
        e0    = 0.5*(nA-2*HtAHt+nHtH);
        
        if strcmp(options.display,'on')==1
            display(logdisplay_parameters(n,r,e0,options));
        end

        % Algorithm
        [H,et,t] = symNMFmex(A,H,maxiter,timelimit,sc);
        
        % Objective function
        % The length of 'et' and 't' (coming from 'symNMFmex')
        % is 'maxiter'. If the timelimit was the first stopping critertion,
        % the unused components are filled with '-1'.
        % Each component et(i) is the decrease of the objective function
        % from et(i-1).
        for i=1:length(et)+1
            if i==1
                e(i) = e0;
            else
                if(et(i-1)~=-1)
                    e(i) = e(i-1)-0.5*et(i-1);
                else
                    break;
                end
            end
        end
        t=[0;t(1:length(e)-1)];
        e=e';
        
        %%%%%%%%%%%% Final objective function %%%%%%%%%%%%%%
        if strcmp(options.display,'on')==1
            nHtH  = norm(H'*H,'fro')^2;
            HtAHt = sum(sum((H'*A).*(H')));
            ef    = 0.5*(nA-2*HtAHt+nHtH);
            display(sprintf('Final objective function=%1.5g',ef));
        end
    end
end

function options=update_options(optnew)
    if nargin==0
        loadoptionsSymNMF;
    else
        loadoptionsSymNMF;
        if isfield(optnew,'maxiter')
            options.maxiter=optnew.maxiter;
        end
        if isfield(optnew,'timelimit')
            options.timelimit=optnew.timelimit;
        end
        if isfield(optnew,'display')
            options.display=optnew.display;
        end
        if isfield(optnew,'shuffle_columns')
            options.shuffle_columns=optnew.shuffle_columns;
        end
        if isfield(optnew,'initmatrix')
            options.initmatrix=optnew.initmatrix;
        end
        if isfield(optnew,'seed')
            options.seed=optnew.seed;
        end
    end
end

function s=logdisplay_parameters(n,r,e0,options)
    s = sprintf('Factorizing a %dx%d matrix using r=%d (maxiter=%d, timelimit=%f)',n,n,r,options.maxiter,options.timelimit);
    switch options.initmatrix
        case 'zeros'
            s = sprintf('%s\nInitial matrix: zeros(%d,%d)',s,n,r);
        case 'dense01'
            s = sprintf('%s\nInitial matrix: rand(%d,%d) with seed=%f',s,n,r,options.seed);
    end
    if options.shuffle_columns
        s=sprintf('%s\nThe columns are shuffled',s);
    end
    s=sprintf('%s\nInitial objective function=%1.5g',s,e0);
end