% Different initialization stragegies for NMF
% 
% W is an m-by-r nonnegative matrix
% H is an r-by-n nonnegative matrix

function [W,H] = initialization(m,n,options,r)
    switch options.rndtype
        case 'rndcube'
            W = rand(m,r);
            H = rand(r,n);
        case 'sparse10'
            [W,H] = initsparse(m,n,r,1,0);
        case 'sparse11'
            [W,H] = initsparse(m,n,r,1,1);
        case 'sparse01'
            [W,H] = initsparse(m,n,r,0,1);
        case 'sparse00'
            [W,H] = initsparse(m,n,r,0,0);
    end
end

% Initialization with a single non-zero in each column or row of W/H
% crW = 1: one per column, crW = 0: one per row
% crH = 1: one per column, crH = 0: one per row
function [W,H] = initsparse(m,n,r,crW,crH)
    if crW==1
        W = matrixsparse(m,r,1);
    elseif crW==0
        W = (matrixsparse(r,m,1))';
    end

    if crH==1
        H = matrixsparse(r,n,1);
    elseif crH==0
        H = (matrixsparse(n,r,1))';
    end
end
% W is a matrix (m x r) with crW ones in each column
function W = matrixsparse(m,r,crW)
    W = zeros(m,r); 
    for i=1:r
        rp=randperm(m);
        W(rp(1:crW),i)=ones(1,crW); 
    end
end