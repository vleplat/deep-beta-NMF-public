% A single MU update of H for the beta-divergence that aims at minimizing 
% 
%               min_{H>=epsilon} D_beta(X,WH). 
% 
% Input: - m-by-n matrix X>=0
%        - m-by-r matrix W>=0
%        - r-by-n matrix H>=0 (initialization)
%        - scalar beta defining the beta-divergence
%        - epsilon: lower bound on the entries of H 
%           (default: machine epsilon). 
% Output: - H is updated using a single MU 
%         - e reports the error *before* the update of H
%
% Code modified from https://sites.google.com/site/nicolasgillis/code of
% the paper N. Gillis, L. T. K. Hien, V. Leplat and V. Y. F. Tan, 
% "Distributionally Robust and Multi-Objective Nonnegative Matrix 
% Factorization", January 2019, http://arxiv.org/abs/1901.10757

function [H,e] = MUbeta(X,W,H,beta,epsilon)

if nargin <= 4
    epsilon = 2^(-52);
end
if nargout == 1
    [N,D] = ND_MUbeta(X,W,H,beta);
else
    [N,D,e] = ND_MUbeta(X,W,H,beta);
end
% Use the gamma power to ensure monotonicity 
if 1 <= beta && beta <= 2
    H = max(epsilon, H .* (N./(D+eps))); 
else
	if beta < 1
        gamma = 1/(2-beta); 
    else
        gamma = 1/(beta-1); 
    end
    H = max(epsilon, H .* ((N./(D+eps)).^gamma) ); 
end