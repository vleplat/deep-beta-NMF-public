% Linear dimensionality reduction: given M (m-by-n) and r < m, 
% compute a projection u (m-by-r), and return 
% Mr = u^T * M (r-by-n) which reduces the dimension of M. 
% 
% *** Description ***
% Given an m-by-n matrix and a dimension r, perform a linear dimensionality
% reduction to dimension r: Mr = u^T*M is r-by-n. 
% 
% [Mr, u] = lindimred(M,r,dimred) 
%
% ****** Input ******
% M          : m-by-n matrix
% r          : dimension of the space to which the data is projected. 
% dimred     : parameter allowing to choose the dimensionality reduction
%              technique used (for m > r): 
%              dimred=0: use SPA (default). 
%              dimred=1: use the truncated SVD, that is, svds(M,r). 
%              dimred=2: use the truncated SVD, but computation based on MM'.
%              else    : use a random projection (u = randn(m,r))
%              (default = 1.)
%
% ****** Output ******
% Mr          : r-by-n reduced data matrix M. 
% u           : corresponding projection, that is, Mr = u'*M. 
% 
% Code adapted from https://sites.google.com/site/nicolasgillis/code

function [Mr, u] = lindimred(M,r,dimred) 

[m,n] = size(M); 
if r > min(m,n)
    error('r has to be smaller than m and n'); 
end
if nargin <= 2, 
    dimred = 0; 
end
if m == r % dimensionality reduction not applicable
    u = eye(r);
else   
    if dimred == 0 % Use SPA
        options.display = 0; 
        K = SPA(M,r,options); 
        u = pinv(M(:,K))'; 
    elseif dimred == 1 % Using truncated SVD
        [u,s,v] = svds(M,r); 
    elseif dimred == 2 % Using truncated SVD of MM'
        MMt = M*M'; 
        [u,s,v] = svds(MMt,r); 
    else  % Using random projection
        u = randn(m,r);  
    end
end
Mr = u'*M;  