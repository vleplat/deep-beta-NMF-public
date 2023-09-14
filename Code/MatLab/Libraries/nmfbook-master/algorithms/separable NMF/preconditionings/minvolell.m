% Active-Set Method for Minimum Volume Ellipsoid
% 
% *** Description ***
% Given the matrix M, compute an approximation of the minimum volume 
% ellipsoid containing the columns of M. 
%
% See Algorithm 2 in N. Gillis and S.A. Vavasis, "Semidefinite Programming 
% Based Preconditioning for More Robust Near-Separable Nonnegative Matrix 
% Factorization", SIAM J. on Optimization 25 (1), pp. 677-698, 2015.
% 
% Q = minvolell(M,delta,addit,affi) 
%
% ****** Input ******
% M          : m-by-n matrix
% delta      : precision 
%              (default = 1e-6, recommended: between 1e-3 and 1e-6.)
% addit      : number of active constraints additional to the r*(r+1)/2
%              (default = r.)
% affi       : display the progess of the algorithm 
%              (default = 1.) 
% 
% ****** Output ******
% Q          : { x | x^T (Q^TQ) x <= 1} is an approximation of the
%              minimum volume ellipsoid containing the columns of M
%              with m_i^T (Q^TQ) m_i <= 1 + delta for all i. 

function Q = minvolell(M,delta,addit,affi) 

[m,n] = size(M); 
r = m; 
if nargin <= 1
    delta = 1e-6;
end
if nargin <= 2
    addit = r; 
end
if nargin <= 3
    affi = 1; 
end

% Number of active columns (could/should be reduced for large r)
Nbr = r*(r+1)/2; 
if Nbr+addit > n
    Nbr = n;
    addit = 0;
end 

etim = cputime; 
if affi == 1
    fprintf('Solving the SDP with precision %0.1d using active set started...\n', delta) 
end

% Inital active set computed using SPA
ind = SPAselect(M,Nbr+addit); 
e = -1; 
k = 1; 
options.display = 0; 
% Main loop
while min(e) <= -delta && k <= 100
    % Create and solve the model 
    cvx_begin 
        cvx_quiet(true)
        variable Q(m,m) symmetric 
        maximize( det_rootn( Q ) )
        subject to
            norms( Q * M(:,ind) , 2 ) <= 1;
    cvx_end
    e = 1-norms( Q * M , 2 ); % Should be nonnegative
    
    % Keep indices with e < 1e-3
    ind = ind(e(ind) <= 1e-3);  
    % If more than r*(r+1)/2 indices: keep only r*(r+1)/2
    if length(ind) > Nbr
       % Using SPA
       indi = SPA(M(:,ind),r,options); 
       % and based on e
       ei = e(ind); ei(indi) = Inf; 
       [a,b] = sort(ei); 
       indi = unique([ indi b( 1:Nbr-r) ]);  
       ind = ind(indi); 
    end

    % Add indices with smallest e (<= 0)
    et = e; et(ind) = 0; 
    for i = 1 : Nbr+addit-length(ind)
        [a,b] = min(et); 
        if a <= 0
            ind = [ind b]; 
            et(b) = 0;
        else
            break;
        end
    end
    
    ind = sort(ind); 
    if affi == 1
        fprintf('%1.0f...',k); 
    end
    k = k+1; 
end

if affi == 1
    fprintf('Solving the SDP required %2.0f changes of active sets. \n',k); 
    fprintf('Solving the original SDP took %2.2f seconds. \n',cputime-etim); 
end

if min(e) <= -delta && k == 101
    warning('The active-set methid did not convege.');
end