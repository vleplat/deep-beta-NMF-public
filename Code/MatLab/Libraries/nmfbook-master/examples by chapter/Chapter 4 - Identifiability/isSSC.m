% This code checks if a matrix satisfies the sufficiently scattered 
% condition (SSC). It only works for small matrices as it relies on vertex 
% enumeration. 
%
% Input : an r-by-n matrix H
% Output: ssc1 = 1 if SSC1 is satisfied, = 0 otherwise.
%         ssc2 = 1 if SSC2 is satisfied, = 0 otherwise.
%            vert is the vertex of H^T x >= 0, e^T x = 1 that maximizes 
%            ||x||^2, that is, vert is an optimal solution of
%            max_x ||x||^2 such that H^T x >= 0, e^T x = 1.
%            If SSCnec=0, vert = [].  
%         SSCnec = 1 if the necessary condition for SSC1 is satisfied, 
%                = 0 otherwise.  
%            SSCnec=1 ensures that {x | H^T x >= 0, e^T x = 1} is bounded.        

function [ssc1,ssc2,vert,SSCnec] = isSSC(H)

% First check the necessary condition for SSC1; see SSC1_nec_cond.m. 
% This condition ensures that { x | H^T x >= 0, e^T x = 1 } is bounded.
SSCnec = SSC1_nec_cond(H);
if SSCnec == 0
    ssc1 = 0; 
    ssc2 = 0; 
    vert = []; 
else
    [r,n] = size(H);
    % Compute the vertices of the dual of cone(H) within e^T x=1.
    V = computevert(H');
    % Look at their norms
    normV = sum(V.^2);
    % Among these vertices, you should have the unit vectors
    ind = find(normV >= 1);
    % Pick vertex that maximizes normV whil having the smallest possible entry
    % (if it is negative, then ssc2 not satisfied)
    [a,b] = min( min( V(:,ind) ) );
    vert = V(:,ind(b));
    if max(normV) > 1+1e-9
        ssc1 = 0;
        ssc2 = 0;
    else
        ssc1 = 1;
        ssc2 = 1; 
        if length(ind) > r
            ssc2 = 0;
        end
    end
end
