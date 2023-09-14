% This codes checks whether a necessary condition for the first
% sufficiently scattered condition (SSC1) is satisfies. 
% 
% More precisely, it checks whether all vectors e-e_k belong to the
% relative intgerior of conv(H); see Algorithm 4.2. and Theorem 4.7 in the 
% book.  
%
% Input : a nonnegative r-by-n matrix H 
% Output: SSCnec = 1 if the necessary condition for SSC1 is satisfied, 
%                = 0 otherwise. 
% 
% Warning: this code can have problem checking the condition if many
% columns of H are very close to one another, making it hard to detect
% numerically whether a point is a vertex or not; this is related to the
% choice of the parameter thres, which is set to 1e-8. 

function SSCnec = SSC1_nec_cond(H) 

[r,n] = size(H); 
if min(H(:)) < 0
    warning('H has negative entries.'); 
    SSCnec = 0;
    return;
end
if rank(H) < r 
    warning('rank(H) < r. The input matrix should be r-by-n with r <= n and rank(H) = r.'); 
    SSCnec = 0;
    return;
end
thres = 1e-6; % Threshold to determine whether a point is a vertex of the 
              % convex hull of given set of points
% Remove zero columns in H 
sumH = sum(H); 
H = H(:,sumH > 0); 
n = size(H,2); 
H = H./repmat(sum(H), r , 1); 
% Remove duplicated columns 
Hn2 = H./repmat(sqrt(sum(H.^2)), r , 1); 
anglesH = triu(Hn2'*Hn2,1); 
[row,col] = find(anglesH >= 1-100*thres); % Remove columns that are too close 
                                          % to one another
ac = unique(col);  
H = H(:, setdiff(1:n,ac)); 
% Check nec condition number of zeros in each row of H is r-1
if min( sum(H'==0) ) < r-1
    SSCnec = 0;
    return;
end
n = size(H,2); 
% Loop to check the condition e-e_k belongs to the relative interior H
for k = 1 : r
    Ik = find(H(k,:) == 0); 
    if rank(H(:,Ik)) < r-1
        SSCnec = 0;
        return;
    else
        % Select vertices
        Ikb = []; 
        for i = 1 : length(Ik)
            C = H(:,Ik([1:i-1 i+1:end])); 
            d = H(:,Ik(i));      
            [~,resnorm] = lsqnonneg(C,d); 
            indi = [1:i-1 i+1:n]; 
            if resnorm > thres
                Ikb = [Ikb Ik(i)]; 
            end
        end
        % check whether e-e_k belong to the relative interior
        d = ones(r,1); d(k) = 0; 
        [x,resnorm] = lsqnonneg(H(:,Ikb),d); 
        if resnorm > thres % || sum( x > thres) < r-1 
        % remove second condition which should only be used with an interior-point method; cf. the book
            SSCnec = 0;
            return;
        end
    end
end
SSCnec = 1; 
