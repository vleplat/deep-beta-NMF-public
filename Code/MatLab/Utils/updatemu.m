function [mu]=updatemu(Phi,Eta,Psi,W,mu,epsi)

[F,K]=size(W);
JF1=ones(F,1);
doLoop = true;
maxitermu = 10^3;
k=1;
while doLoop && k<=maxitermu
    mu_prev=mu;
    Mat=W .*(((Phi+JF1*mu').^2+Eta).^(1/2)-(Phi+JF1*mu'))./(Psi+eps);
    xi=(sum(Mat,1)-ones(1,K))';
    Matp=(W./Psi+eps).*((Phi+JF1*mu')./(((Phi+JF1*mu').^2+Eta).^(1/2))-ones(F,K));
    xip=sum(Matp,1)';
    mu=mu-xi./xip;
    if(max(abs(mu-mu_prev))<=epsi)
        doLoop=false;
    end
    k=k+1;
end

if k == maxitermu
    warning('Newton-Raphson procedure reached maximum of iterations.');
end
% flag=1; %uncomment for debugging and check the convergence rate of mu
% figure;
% semilogy(max(xi_save,0)')

end%EOF