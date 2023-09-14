% Find vertices V of the set P = { x in R^{k-1} | Cx+d >= 0 } with brute force
function V = vertices(C,d);

[m,k] = size(C);
V = []; lP = 0;
choices = nchoosek(1:length(C(:,1)),k);
% Choose k lin. ind. inequalities from Cx+d >= 0 and compute the intersection
for i = 1 : length(choices(:,1))
    if rank(C(choices(i,:),:)) == k
        x = C(choices(i,:),:)\[-d(choices(i,:))];
        % Check if the intersection is in P
        if  min(C*x+d) >= -1e-9 && (lP == 0 || min(sum((V-repmat(x,1,lP)).^2))>1e-6)
            V = [V x]; lP = lP+1;
        end
    end
end