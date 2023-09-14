% Generates the following linear Euclidean matrix of dimension n: 
% M(i,j) = (i-j)^2. 
% See `Heuristics for Exact Nonnegative Matrix Factorization', A. Vandaele,  
% N. Gillis, F. Glineur and D. Tuyttens, for more details. 

function [S,r]=LEDM(n)

    S = (repmat(1:n,n,1)-repmat((1:n)',1,n)).^2;
    
    % For n <= 5, rank_+ = n. 
    % Gillis, Glineur, `On the Geometric Interpretation of the Nonnegative 
    % Rank', Linear Algebra and its Applications 437(11):2685-2712, 2012.
    % 
    % Moreover, we have the following upper bounds
    %     rank_+(n) <= rank_+(n/2) + 2 for n even, and 
    %     rank_+(n) <= rank_+(n+1). 
    % See P. Hrubes: On the nonnegative rank of distance matrices, 
    % Information Processing Letters 112(11):457�461, 2012. 
    % 
    % We conjecture the the upper bound obtained combining the results in 
    % the two papers above is tight: 
    nn = n; k = 0; 
    while nn > 5
        nn = nn/2;
        k = k + 1;
    end
    r = ceil(nn) + 2*k; 
end