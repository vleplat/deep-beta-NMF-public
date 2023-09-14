% Normalization the pair (W,H) depending on the NMF model. 
% Four possibilites: 
% sumtoone = 1: scale H such that H^Te <= e, and update W accordingly. 
% sumtoone = 2: scale H such that He = e, and update W accordingly. 
% sumtoone = 3: scale H such that W^Te = e, and update H accordingly. 
% sumtoone = 4: scale H such that H^Te = e, and update W accordingly. 
% 
% Note: in the cases sumtoone in {1,4}, this is not wlog and could modify 
% the product W*H. (In fact, we can only modify the row sums of H or the
% column sums of W wlog, this is the scaling degree of freedom in NMF.) 

function [W,H] = normalizeWH(W,H,sumtoone,X) 

if sumtoone == 1 % Normalize so that H^Te <= e
                 % entries in cols of H sum to at most 1
    Hn = SimplexProj( H );
    if norm(Hn - H) > 1e-3*norm(Hn); 
       H = Hn; 
       % reoptimize W, because this normalization is NOT w.l.o.g. 
       options.inneriter = 100; 
       options.H = W'; 
       W = nnls_FPGM(X',H',options); 
       W = W'; 
    end
    H = Hn; 
elseif sumtoone == 2 % Normalize so that He = e, 
                     % entries in rows of H sum to 1
    scalH = sum(H');
    H = diag( scalH.^(-1) )*H;
    W = W*diag( scalH );
elseif sumtoone == 3 % Normalize so that W^T e = e, 
                     % entries in cols of W sum to 1
    scalW = sum(W);
    H = diag( scalW )*H;
    W = W*diag( scalW.^(-1) );
elseif sumtoone == 4 % Normalize so that H e = e, 
                     % entries in cols of H sum to 1  
    Hn = SimplexColProj( H );
    if norm(Hn - H) > 1e-3*norm(Hn); 
       H = Hn; 
       % reoptimize W, because this normalization is NOT w.l.o.g. 
       options.inneriter = 100; 
       options.H = W'; 
       W = nnls_FPGM(X',H',options); 
       W = W'; 
    end
    H = Hn;   
end