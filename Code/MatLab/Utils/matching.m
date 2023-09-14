function [W1,H1,W2,H2]=matching(W1,H1,W2,H2,time)
[~, K]=size(W1);

% Normlization of rows of Hi
D1=diag(1./sum(H1,2));
H1_norm=D1*H1;

D2=diag(1./sum(H2,2));
H2_norm=D2*H2;

% Creation of the cost matrix Cij=|H1(i,:)-H2(j,:)|_1
C=zeros(K,K);
for i=1:K
    for j=1:K
        C(i,j)=norm(H1_norm(i,:)-H2_norm(j,:),1);
    end
end

[assignment,~] = munkres(C);

% Creation of the permutation Matrix
P=eye(K,K);
I=eye(K,K);
for k=1:K
    P(k,:)=I(assignment(k),:);
end

H2=P*H2;
W2=W2*P';

end%EOF