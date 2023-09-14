% Heuristic Multi-Start 1 for Non-Negative Matrix Factorization
% 
% [W,H,e,e0] = heuristic_ms1(X,r,options)
%
% Input.
%   X       : (m x n) matrix to factorize
%   r       : the non-negative rank
%   options : set of parameters (options is required).
% Output.
%   (W,H)   : nonnegative matrices s.t. WH approximates X
%   e       : error of the solution (W,H)
%   e0      : error of the solution before Local Improvement

function [W,H,e,e0] = heuristic_ms1(X,r,options)
    [W,H] = initialization(size(X,1),size(X,2),options,r);
    e0    = norm(X-W*H,'fro')/norm(X,'fro');
    if e0>options.tolerance
        [W,H,e] = localimprovement(X,W,H,options);
    else
        e=e0;
    end
end