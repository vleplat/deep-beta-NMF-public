% Generating synthetic data sets - 4 types
%
% W is m-by-r
% H = [I_r, H'] is r-by-(n+r)
% N is m-by-(n+r)
% delta is the noise level
% xp is the type of data set:
%   =1: well-conditioned Dirichlet
%   =2: well-conditioned Middle points
%   =3: ill-conditioned Dirichlet
%   =4: ill-conditioned Middle points
% for xp=3,4, cond(W) = 10^condW

function [W,H,Noise] = synthdatasetSepNMF(m,n,r,delta,xp,condW,diri)

if nargin <= 6
    diri = 0.5;
end
if nargin <= 5
    condW = 3;
end
if r > m
    warning(' The ambiant dimension m is smaller than r.');
end
% Generate W
gamma2 = 0;
while gamma2 < 10^(-condW) % This loop is useful when r < m, 
                           % In fact, gamma2(W) > sigma_r(W)
    W = rand(m,r);
    if xp >= 3 % ill-conditioned
        [u,s,v] = svds(W,r);
        W = u*diag( logspace(-condW,0,r) )*v';
    end
    gamma2 = gamma2param(W);
end
if nargout >= 2
    % Generate H and N
    % Dirichlet experiment
    if xp == 1 || xp == 3
        % H is generated so that the columns of W are repeated twice
        % and the remaining columns are drawn following a Dirichlet
        % distribution
        if diri > 0
            alpha = diri*ones(r,1);  % Parameter of the Dirichlet distribution
        else
            alpha = rand(r,1);  % Parameter of the Dirichlet distribution
        end
        H = [eye(r) sample_dirichlet(alpha,n)'];
        Noise = [randn(m,n+r)];  % Noise is Guassian
        Noise = delta * Noise/norm(Noise,'fro') * norm(W*H,'fro');
        % Middle points experiment
    elseif xp == 2 || xp == 4
        % H is generated so that the first r columns of M are the columns
        % of W, and the remaining ones are on the middle point of any
        % combination of two columns of W
        n = nchoosek(r,2);
        H = [eye(r) nchoose2(r)/2];
        X = W*H;
        Noise = delta*[zeros(m,r) X(:,r+1:end)-repmat(mean(W,2),1,n)];
    end
end