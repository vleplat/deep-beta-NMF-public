% symNMF on Karate club; see Figure 5.9
% 
% To display the results, we modified the code available from 
% https://blogs.mathworks.com/loren/2015/09/30/can-we-predict-a-breakup-social-network-analysis-with-matlab/

clear all; close all; clc;

load('karate.mat');
% Construct X 
edg = size(edges,1); 
n = max(edges(:)); 
X = eye(n); 
for i = 1 : edg 
    X(edges(i,1), edges(i,2)) = 1;
    X(edges(i,2), edges(i,1)) = 1;
end

% Inner rank of the factorization
r = 2;

% Options of symNMF (see also loadoptions.m)
options.maxiter   = 1000;
options.timelimit = 5;
options.initmatrix='dense01'; % random initialization
options.seed=0; % Use same seed to have the same initialization 
% Running symNMF 
options.shuffle_columns = 0;
[W,e,t] = symNMF(X,r,options); % The model is X ~ WW^T. 
% Display result 
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 2);
figure; 
subplot(1,3,1);
imagesc(X); 
title('X'); 
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
subplot(1,3,2);
imagesc(W(:,1)*W(:,1)'); 
title('W(:,1)*W(:,1)^T'); 
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
subplot(1,3,3);
imagesc(W(:,2)*W(:,2)'); 
title('W(:,2)*W(:,2)^T'); 
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
W 
% Plot graphs 
intes = 4; 
X = X-eye(n); 
G = graph(X);           % create a graph from edges
G.Nodes = table(name);  % name the nodes
figure; 
subplot(1,3,1);  
plot(G,'LineWidth',intes*G.Edges.Weight);
title('X'); 
% Add weights from H(:,1)*H(:,1)' 
WWt1 = W(:,1)*W(:,1)'+1e-6; 
for i = 1 : edg 
    G.Edges.Weight(i) = WWt1( G.Edges.EndNodes(i,1), G.Edges.EndNodes(i,2) );
end
subplot(1,3,2); 
p1 = plot(G,'LineWidth',intes*G.Edges.Weight); 
title('W(:,1)*W(:,1)^T'); 
for i = 1 : edg
    if G.Edges.Weight(i) < 0.1 
        highlight(p1,[G.Edges.EndNodes(i,1), G.Edges.EndNodes(i,2)],'EdgeColor','r')
    end
end
% Add weights from H(:,2)*H(:,2)' 
WWt2 = W(:,2)*W(:,2)'+1e-6; 
for i = 1 : edg 
    G.Edges.Weight(i) = WWt2( G.Edges.EndNodes(i,1), G.Edges.EndNodes(i,2) );
end
subplot(1,3,3);  
p1 = plot(G,'LineWidth',intes*G.Edges.Weight); 
title('W(:,2)*W(:,2)^T'); 
for i = 1 : edg
    if G.Edges.Weight(i) < 0.1 
        highlight(p1,[G.Edges.EndNodes(i,1), G.Edges.EndNodes(i,2)],'EdgeColor','r')
    end
end