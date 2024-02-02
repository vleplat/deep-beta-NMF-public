function [mu]=updatemu_hcols(C,D,H,gamma_beta,mu,epsi)
[K,T]=size(C);
JK1=ones(K,1);
delta=1;
doLoop = true;
maxitermu = 10^4;
k=1;
while doLoop && k<=maxitermu
    mu_prev=mu;
    % Mat=H .* (C./(D-JK1*mu'+eps)).^(gamma_beta);
    Mat=H .* (C./(D-JK1*mu'+eps));
    xi=(sum(Mat,1)-delta*ones(1,T))';
    % Matp=H .* (C./(D-JK1*mu'+10^-8)).^(gamma_beta-1);
    Matp=H;
    Matp=1*Matp.*C./(D-JK1*mu'+eps).^2;
    xip=sum(Matp,1)';
    mu=mu-xi./(xip);
    if(max(abs(mu-mu_prev))<=epsi)
        doLoop=false;
    end
    k=k+1;
end

end%EOF