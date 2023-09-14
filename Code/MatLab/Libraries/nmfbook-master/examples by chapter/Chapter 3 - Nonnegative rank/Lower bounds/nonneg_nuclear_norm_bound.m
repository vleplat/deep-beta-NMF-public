% Lower bound for the nonnegative rank using the nonnegative nuclear norm
% approach proposed in 
% Fawzi, H., & Parrilo, P. A. (2015). Lower bounds on nonnegative rank via 
% nonnegative nuclear norms. Mathematical Programming, 153(1), 41-66.
% 
% This code uses CVX that you need to download from http://cvxr.com/cvx/
% 
% Input : matrix X 
% Output: lower bound lb <= rank_+(X) 
%         and solution (nu,Y,Z) of the corresponding semidefinite program
% See the paper for more details, or Section 3.4.8 of the book. 

function [lb,nu,Y,Z] = nonneg_nuclear_norm_bound(X) 

[m,n] = size(X); 

cvx_begin sdp quiet
    variable Y(m,m);
    variable Z(n,n);
    variable B(m+n,m+n);
          
    minimize( 0.5* (trace(Y) + trace(Z))  );
    subject to
        B == [ Y X ; X' Z ]; 
        B == hermitian_semidefinite( m+n ); 
        B(:) >= 0; 
cvx_end

nu = 0.5* (trace(Y) + trace(Z)); 
lb = (nu/norm(X,'fro'))^2; 