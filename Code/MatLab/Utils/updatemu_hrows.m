function [mu]=updatemu_hrows(C,D,H,gamma_beta,mu,epsi)
[K,T]=size(C);
% JK1=ones(K,1);
JN1 = ones(T,1);
delta=1;
doLoop = true;
maxitermu = 10^4;
k=1;
while doLoop && k<=maxitermu
    mu_prev=mu;
    Mat=H .* (C./(D-mu*JN1'+10^-8)).^(gamma_beta);
    xi=(sum(Mat,2)-delta*ones(K,1));
    Matp=H .* (C./(D-mu*JN1'+10^-8)).^(gamma_beta-1);
    Matp=1*Matp.*C./(D-mu*JN1'+10^-8).^2;
    xip=sum(Matp,2);
    mu=mu-xi./(xip);
    if(max(abs(mu-mu_prev))<=epsi)
        doLoop=false;
    end
    k=k+1;
end

if k == maxitermu
    warning('Newton-Raphson procedure reached maximum of iterations.');
end

end%EOF