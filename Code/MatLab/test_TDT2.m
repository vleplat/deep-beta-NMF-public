% test on TDT2 data set 
clear all; clc; 
load tdt2_top30 
X = X'; 
rng(2020); 
r = 20; 
[W,H] = FroNMF(X,r);
for i = 1 : r
    [a,b] = sort(W(:,i),'descend');
    for j = 1 : 10 % Keep the 10 words with the largest value in W(:,i) 
        Topics{j,i} = words{b(j)};
    end
end 
% Subampled words (rows of X)  
K = []; 
maxowrdspertop = 500/r; 
for i = 1 : r
    [~,coli] = sort(W(:,i),'descend'); 
    K = unique([K; coli(1:maxowrdspertop)]); 
end

Xk = X(K,:); 
wordsk = words(K); 
[Wk,Hk] = FroNMF(Xk,r);
for i = 1 : r
    [a,b] = sort(Wk(:,i),'descend');
    for j = 1 : 10 % Keep the 10 words with the largest value in W(:,i) 
        Topicsk{j,i} = wordsk{b(j)};
    end
end 

% %% Deep NMF 
% r = [20 10]; 
% options.delta = ones(1,length(r));
% options.alpha_tilde = 0.001*[1; 1];
% 
% [Wl,Hl,el,inWH,output] = deepKL_NMF(Xt',r,options);
% 
% for i = 1 : r(1)
%     [a,b] = sort(Wl{1}(:,i),'descend');
%     for j = 1 : 10 % Keep the 10 words with the largest value in W(:,i) 
%         Topics{j,i} = wordsK(b(j));
%     end
% end 