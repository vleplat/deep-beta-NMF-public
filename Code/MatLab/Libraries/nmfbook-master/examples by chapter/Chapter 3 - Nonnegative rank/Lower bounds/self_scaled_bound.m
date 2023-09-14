% Lower bound for the nonnegative rank using the self-scaled bound
% proposed in 
% Fawzi, H., & Parrilo, P. A. (2016). Self-scaled bounds for atomic cone 
% ranks: applications to nonnegative rank and cp-rank. Mathematical 
% Programming, 158(1-2), 417-465.
% 
% This code uses CVX that you need to download from http://cvxr.com/cvx/
% 
% Input : matrix X 
% Output: lower bound tausos <= rank_+(X) 
%         and solution Y of the corresponding semidefinite program
% See the paper for more details, or Section 3.4.9 of the book. 

function [tausos,Y] = self_scaled_bound(X) 

[m,n] = size(X); 

cvx_begin sdp quiet
    variable Y(m*n,m*n);
    variable tausos(1);
          
    minimize( tausos );
    subject to
        [ tausos vec(X)' ; vec(X) Y ] == hermitian_semidefinite( m*n+1 ); 
        for i = 1 : m
            for j = 1 : n
                Y(i+(j-1)*m,i+(j-1)*m) <= X(i,j)^2; 
                for k = i+1 : m
                    for l = j+1 : n
                        Y(i+(j-1)*m,k+(l-1)*m) == Y(i+(l-1)*m,k+(j-1)*m); 
                    end
                end
            end
        end
cvx_end

