% Post-processed successive projection algorithm (Post-SPA) which is 
% equivalent to FastAnchorWords under the assumptions that the sum of the 
% entries of each columns of H is smaller than one instead of being equal 
% to one. 
% 
% *** Description ***
% At each step of the algorithm, the column of M maximizing ||.||_2 is 
% extracted, and M is updated by projecting its columns onto the orthogonal 
% complement of the extracted column (see FastSepNMF.m). 
% After r indices have been identified, the post-processing of Arora et al. 
% (A Practical Algorithm for Topic Modeling with Provable Guarantees, 
% ICML '13) is used to refine the solution. 
% 
% J = FastAnchorWords(M,r,normalize) 
%
% ****** Input ******
% M = WH + N : a (normalized) noisy separable matrix, that is, W full rank, 
%              H = [I,H']P where I is the identity matrix, H'>= 0 and its 
%              columns sum to at most one, P is a permutation matrix, and
%              N is sufficiently small. 
% r          : number of columns to be extracted. 
% Options: 
% .normalize  : normalize=1 will scale the columns of M so that they sum to one,
%              hence matrix H will satisfy the assumption above for any
%              nonnegative separable matrix M. 
%              normalize=0 is the default value for which no scaling is
%              performed. For example, in hyperspectral imaging, this 
%              assumption is already satisfied and normalization is not
%              necessary. 
% .display    : =1 displays the prgression of the iterations
%
% ****** Output ******
% J        : index set of the extracted columns. 

function J = FastAnchorWords(M,r,options) 

[m,n] = size(M); 
% Options 
if nargin <= 2
    options = [];
end
if ~isfield(options,'normalize')
    options.normalize = 0; 
end
if ~isfield(options,'display')
    options.display = 1; 
end

if options.normalize == 1
    % Normalization of the columns of M so that they sum to one
    D = spdiags((sum(M).^(-1))', 0, n, n); 
    M = M*D; 
end
if options.display == 1
    fprintf('Extraction of the indices by FastAnchorWords: \n'); 
end
% First phase = Successive Projection Algorithm
optionsSPA.display = 0; 
J = SPA(M,r,optionsSPA);   

% Second phase = Post-Processing of Arora et al. 
r = length(J); 
for j = 1 : r
    [k,normM,U] = SPArefine(M,J([1:j-1 j+1:r])); 
    J(j) = k; 
    if options.display == 1
        fprintf('%2.0f...', j); 
        if mod(j,10) == 0
            fprintf('\n');
        end
    end
end
if options.display == 1
    fprintf('\n'); 
end