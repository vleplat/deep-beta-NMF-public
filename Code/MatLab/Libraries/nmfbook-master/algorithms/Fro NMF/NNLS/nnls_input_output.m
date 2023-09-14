% Input/output for NNLS algorithms, nnls_x  
% 
% ****** Input ******
%   X  : m-by-n matrix
%   W  : m-by-r matrix
% 
% ---Options--- 
%   .init        : r-by-n initialization matrix. 
%                   default = [] leading to the use of nnls_init.m. 
% .delta in [0,1): Stop iterations when the distance between two iterates
%                 is smaller than the first two ones * delta. 
%                   defaut = 1e-6.
% .inneriter >= 1: maximum number of inner itrations. 
%                   default = 500.
%
% ****** Output ******
%   H  : an r-by-n nonnegative matrix which is approximate solution of 
%           min_{H >= 0} ||X-WH||_F^2
% WTW  : W'*W, and  
% WTX  : W'*X which can be used to compute the error efficiently since 
%        ||X - WH||_F^2 = ||X||_F^2 - 2 <WTX,H> + <WTW,H*H^T>