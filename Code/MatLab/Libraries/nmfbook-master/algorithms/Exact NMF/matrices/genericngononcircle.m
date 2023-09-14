% Generates the slack matrix of a generic n-gon with vertices
% lying on a circle. The matrix S has dimension (n x n).
% See `Heuristics for Exact Nonnegative Matrix Factorization', A. Vandaele,  
% N. Gillis, F. Glineur and D. Tuyttens, for more details.

function S = genericngononcircle(n)
%     d = 3;
%     temp = 2*pi/(n*d)*(2*rand - 1);
%     if temp<0
%         temp = 2*pi+temp;
%     end
%     angles(1,1) = temp;
%     for i=1:n-1
%         a = 2*pi*i/n - 2*pi/(d*n);
%         b = 2*pi*i/n + 2*pi/(d*n);
%         angles(i+1,1) = a + (b-a)*rand;
%     end

    angles = rand(n,1)*2*pi; 
    P = [cos(angles) sin(angles)];
    
    o=convhull(P,'simplify',true);
    V=P;%-ones(n,1)*mean(P);
    plot(P(:,1),P(:,2),'red+'); hold on;
    plot(V(o,1),V(o,2),'o-'); figure(gcf);
    for i=1:n
        F(i,1:2) = [V(o(i),2)-V(o(i+1),2) V(o(i+1),1)*V(o(i),1)] / (V(o(i),2)*V(o(i+1),1)-V(o(i),1)*V(o(i+1),2));
        v1=V(o(i),1); v2=V(o(i),2);
        w1=V(o(i+1),1); w2=V(o(i+1),2);
        F(i,1:2) = [v2-w2 w1-v1]/(v2*w1-v1*w2);
        %buggy: H = [V(o(i),2)-V(o(i+1),2) V(o(i+1),1)*V(o(i),1)] / (V(o(i),2)*V(o(i+1),1)-V(o(i),1)*V(o(i+1),2));
    end
    S=max(ones(n)-F*V',0);
end