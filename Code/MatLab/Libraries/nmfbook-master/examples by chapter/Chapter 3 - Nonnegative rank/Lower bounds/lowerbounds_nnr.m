% Compare lower bounds for the nonnegative rank on a m-by-n matrix X
%
% Input : an m-by-n nonnegative matrix X, with rank_+^*(X) = rnrank.
% Output: lower bounds for the nonnegative rank of X:
%         - rc      : rectangle covering bound of X
%         - geo     : geometric bound, requires the knowledge of rnrank
%                     (if it is not provided, it is not computed)
%         - nnucnorm: nonnegative nuclear norm bound
%         - tausos  : sum-of-squares bounds
%         - hypse   : hyperplane separation bound

function [rc,geo,nnucnorm,tausos,hypsep] = lowerbounds_nnr(X,rnrank)

% Rectangle covering bound
disp('*****************************************************************************');
disp(' Rectangle covering bound via linear integer programming');
disp('*****************************************************************************');
options = optimoptions('intlinprog','Display','off');
rc = rec_cov_bound(X,options)

fprintf('Click any button to continue.');
pause;
fprintf('\n');

if nargin >= 2
    % Geometric bound, requires rnrank as an input 
    disp('*****************************************************************************');
    disp(' Lower bound based on the geometric interpretation (Gillis & Glineur, 2012)');
    disp('*****************************************************************************');
    geo = geometric_bound(rank(X),size(X,2))
    
    fprintf('Click any button to continue.');
    pause;
    fprintf('\n');
else
    geo = 0; 
end

% Nuclear norm bound
disp('!!! For the next bound, you need to install CVX !!!');
disp('!!! See http://cvxr.com/cvx/                    !!!');
disp('*****************************************************************************');
disp(' Lower bound based on the nonnegative nuclear norm (Fawzi & Parrilo, 2015)');
disp('*****************************************************************************');
nnucnorm = nonneg_nuclear_norm_bound(X)

fprintf('Click any button to continue.');
pause;
fprintf('\n');

% Self-scaled bound
disp('*****************************************************************************');
disp(' Lower bound based on the self-scaled bound (Fawzi & Parrilo, 2016)');
disp('*****************************************************************************');
tausos = self_scaled_bound(X)

fprintf('Click any button to continue.');
pause;
fprintf('\n');

% Hyperplane bound
disp('*****************************************************************************');
disp(' Lower bound based on the hyperplane separation bound (Rothvoss, 2015) using');
disp('*****************************************************************************');
Z = X - 0.5 + 1e-6; % This is an arbitrary choice for Z which seems to work 
                    % reasonnably well for the integeer matrices presented 
                    % in the book
indineg = find(Z < 0);
Ztrue = Z;
Z(indineg) = -1000;
Z
[alphaZ,xopt,yopt] = hyperplane_separation_bound(Z);
hypsep = sum(sum( Z.*X ) )/alphaZ/max(X(:))
disp('*****************************************************************************');