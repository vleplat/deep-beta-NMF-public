%      Deep Nonnegative Matrix Factorization with Beta Divergences        %
%-------------------------------------------------------------------------%

% Copyright (c) 2023 Valentin Leplat, Nicolas Gillis, Akwum Onwunta and Le Thi Khanh Hien
% Contact: v.leplat@skoltech.ru or Nicolas.GILLIS@umons.ac.be

% This software reproduces the results from the preprint called:
% "Deep Nonnegative Matrix Factorization with Beta Divergences" - V.
% Leplat, N.Gillis, A. Onwunta and L.-T.-K. Hien

%In order to run the demos, you will need to add to your MATLAB path:

% - NMF toolboxes from Nicolas Gillis book: https://gitlab.com/ngillis/nmfbook

%-------------------------------------------------------------------------%
%                              CONTENT
% - /Libraries : contains helpful libaries.
% - /demos : contains demo files that produce tables and figures.
% - /Datasets_TopicModel : contains data sets for tests reported in Section 4.2 of the preprint
%  - main.m : codes that allow to run desired demos.
%  - /Utils : contains helpful files and MatLab routines to run the demos.
%-------------------------------------------------------------------------%
%                               MENU
% You can launch a specified demo by typing its number. The resulting tables
% and figures produced will be stored in the figures folder.
%
% 1:  Demo for CBCL produces Fig. 2-3 and Table 2
% 2:  Demo for topic modeling, see section 4.2
% 3:  Benchmarking for HSI tests produces Fig. 5-10 
%-------------------------------------------------------------------------%
addpath(genpath('./'));
list_demos = ["test_deepKLNMF_CBCL","test_deepKLNMF_text","test_deepKLNMF_HSI"]; 

disp('Which file do you want to run ?'); 
disp('Enter 1 for facial features extraction (CBCL)'); 
disp('Enter 2 for topic modeling (TDT2-top30)'); 
disp('Enter 3 for hyperspectral unmixing'); 
prompt = "";

num = input(prompt);
eval(list_demos(num));
