% Normalization of the columns of W and rows of H such that 
% ||W(:,k)||_2 = ||H(k,:)||_2 for 1 <= k <= r.

function [W,H]=normNMF(W,H)
    [~,r]=size(W);
    for k=1:r
        alpha_wk = norm(W(:,k));
        alpha_hk = norm(H(k,:));
        W(:,k)   = W(:,k)*sqrt(alpha_wk*alpha_hk)/alpha_wk;
        H(k,:)   = H(k,:)*sqrt(alpha_wk*alpha_hk)/alpha_hk;
    end    
end