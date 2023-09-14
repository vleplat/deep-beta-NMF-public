% Sparse NMF using fast gradient and coordinate descent 
% 
% min_{W>=0,H>=0} || X - WH ||_F such that the average Hoyer sparsity of 
% the columns of W is given by s. 
%
% See the paper Grouped sparse projection by Nicolas Gillis, Riyasat Ohib, 
% Sergey Plis and Vamsi Potluru, http://arxiv.org/abs/1912.03896, 2019. 
% 
% ******
% Input  
% ******
% X      : m-byn matrix to factorize
% r      : rand of the factorization 
% s      : average sparsity of the columns of W, the first factor in the
%          factorization. For s=[], no sparsity constraint enforced. 
% options: 
%          - sW       : average sparsity of the columns of W, the first 
%                       factor in the factorization. For sW=[] (default), no 
%                       sparsity constraint enforced. 
%          - sH       : average sparsity of the rows of H, the second  
%                       factor in the factorization. For sH=[] (default), no 
%                       sparsity constraint enforced. 
%          - maxiter  : maximum number of iterations (number of times W and H
%                       are updated). -default=500
%          - timemax  : maximum alloted time 
%          - delta    : stopping criterion for inner iterations  -default=0.1
%          - inneriter: maximum number of inner iterations when updating W
%                       and H. -default=10
%          - W and H  : initial matrices. default-rand(m,r) and rand(r,n) 
%                       +scaling
%          - w        : weights w{i} i=1,2,...,r to compute the projection
%                       -default=/, meaning w{i}=ones(m,1) for all i
%                       (standard Hoyer sparsity). 
%          - FPGM     : Update of H using FPGM if options.FPGM = 1,
%                       (default:0 and HALS is used)
%          - display  : = 1, displays evolution of the error 
%
% ******
% Output 
% ******
% W (m by r) and H (r by n) such that ||WH-W||_F is small and the averge
% sparsity of the columns of W is s. 
%
% (e,t): error and time, plot(t,e) plots the error ||X-WH||_F over time

function [W,H,e,t] = sparseNMF(X,r,options)

timect = cputime; 
[m,n] = size(X); 
if nargin <= 2
    options = [];
end
if ~isfield(options,'sW')
    options.sW = []; 
end
if ~isfield(options,'sH')
    options.sH = []; 
end
if ~isfield(options,'maxiter')
    options.maxiter = 500; 
end
if ~isfield(options,'timemax')
    options.timemax = 5; 
end
if ~isfield(options,'delta')
    options.delta = 0.1; % Stop inner iter if improvment of current 
                         % iteration is not as good as delta * improvment 
                         % of the first one. 
end
if ~isfield(options,'inneriter')
    options.inneriter = 10; % Maximum number of inner iterations 
end
if ~isfield(options,'FPGM')
    options.FPGM = 0;   % Use of FGM for update of H 
end
if ~isfield(options,'colproj')
    options.colproj = 0;   % Use of column-wise projection for 
end
if isfield(options,'W') && isfield(options,'H') % Initialization 
    W = options.W;
    H = options.H; 
else
    W = rand(m,r); 
    H = rand(r,n); 
    % Optimal scaling  
    alpha = sum(sum( W'*X .* H) ) /  sum(sum( (W'*W).*(H*H') ) ); 
    W = alpha*W; 
end
if ~isfield(options,'display')
    options.display = 1; 
end
% Projection of W to achieve average sparsity s 
if ~isempty(options.sW)
    W = weightedgroupedsparseproj_col(W,options.sW,options);
end
if ~isempty(options.sH)
    H = weightedgroupedsparseproj_col(H',options.sH,options);
    H = H'; 
end
itercount = 1; 
nX2 = sum(X(:).^2); 
nX = sqrt(nX2); 
e = []; 
t = []; 
Wbest = W; 
Hbest = H; 
ndisplay = 10; 
if options.display == 1
    disp('Iteration number and relative error in percent:') 
end
while itercount <= options.maxiter && cputime-timect <= options.timemax
    % Normalize W and H so that columns/rows have the same norm, that is, 
    % ||W(:,k)|| = ||H(k,:)|| for all k 
    normW = sqrt((sum(W.^2)))+1e-16; 
    normH = sqrt((sum(H'.^2)))+1e-16; 
    Wo = W; Ho = H; 
    for k = 1 : r
        W(:,k) = W(:,k)/sqrt(normW(k))*sqrt(normH(k)); 
        H(k,:) = H(k,:)/sqrt(normH(k))*sqrt(normW(k)); 
    end
    % Update H: 
    % (1) No sparsity constraints: block coordinate descent method; 
    % See N. Gillis and F. Glineur, "Accelerated Multiplicative Updates 
    % and Hierarchical ALS Algorithms for Nonnegative Matrix Factorization", 
    % Neural Computation 24 (4), pp. 1085-1105, 2012.
    if options.FPGM == 0 && isempty(options.sH)
        options.init = H; 
        H = NNLS(W,X,options); 
    else
    % (2) With sparsity constraints: fast projected gradient method (FGM) 
    % similar to NeNMF from 
    % Guan, N., Tao, D., Luo, Z., & Yuan, B. NeNMF: An optimal 
    % gradient method for nonnegative matrix factorization, IEEE 
    % Transactions on Signal Processing, 60(6), 2882-2898, 2012.
        options.s = options.sH; 
        H = fastgradsparseNNLS(X',W',H',options); 
        H = H'; 
    end
    % Update W: same as for W
    if isempty(options.sW) % A-HALS
        options.init = W'; 
        [W,HHt,XHt] = NNLS(H',X',options);         
        W = W'; 
        XHt = XHt'; 
    else                   % FPGM 
        options.s = options.sW; 
        [W,XHt,HHt] = fastgradsparseNNLS(X,H,W,options); 
    end
    % Time and error
    e = [e sqrt( max(0, (nX2-2*sum(sum(W.*XHt))+ sum(sum(HHt.*(W'*W)))) ) )/nX]; 
    t = [t cputime-timect]; 
    % Keep best iterate in memory as FGM is not guaranteed to be monotone
    if itercount >= 2 
        if e(end) <= e(end-1)
            Wbest = W; 
            Hbest = H; 
        end
    end
    % Display
    if options.display == 1
        if mod(itercount,ndisplay) == 0
            fprintf('%2.0f:%2.3f - ',itercount,100*e(itercount));
        end
        if mod(itercount,ndisplay*10) == 0
            fprintf('\n');
        end
    end
    itercount = itercount+1; 
end
if mod(itercount,ndisplay*10) > 0, fprintf('\n'); end 
W = Wbest; 
H = Hbest; 