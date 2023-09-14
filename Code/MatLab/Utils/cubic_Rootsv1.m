function [y]=cubic_Rootsv1(a_tilde,b_tilde,c_tilde,d_tilde,y0,maxiter)

%Check that a_tilde is different from 0

if(a_tilde==0)
    %disp('quadratic model')
    del=c_tilde^2-4*b_tilde*d_tilde;
    if(del>=0)
        x1=(-c_tilde+sqrt(del))/(2*b_tilde);
        x2=(-c_tilde-sqrt(del))/(2*b_tilde);
        y=max(x1,x2);
    else
        y=y0;
        %disp('no real roots')
    end
end

%Check if there are three real roots
a=b_tilde/a_tilde;
b=c_tilde/a_tilde;
c=d_tilde/a_tilde;

p=b-3*(a/3)^2;
q=c-a/3*(p+(a/3)^2);
q_m=2*(sign(p)*p/3)^(3/2);
alpha=q/q_m;
if(p<=0 && abs(q)<=q_m)
    %disp('three real roots')
    xm=sqrt(-1/3*p);
    if(abs(alpha)>0 && abs(alpha)<=0.35)
        k=2/3*alpha;
        h_0=k*(1+1/3*k^2+1/3*k^4+4/9*k^6+55/81*k^8);
        phi_alpha=h_0;   %estimate of r2 (intermediate root)
        % Adusting tje initial estimate: Ostrowksy
        u_0=phi_alpha;
        for iter=1:maxiter
            K = Ostrowksy(u_0,alpha);
            u1=u_0-K;
            u_0=u1;
        end
        r2=u1;
        r3=-r2/2+sqrt(12-3*r2^2)/2;
        r1=-r2/2-sqrt(12-3*r2^2)/2;
        
        %Obtaining the roots of the original cubic
        x1=r1*xm;
        x2=r2*xm;
        x3=r3*xm;
        
        y1=x1-a/3;
        y2=x2-a/3;
        y3=x3-a/3;
        
        y=[y1 y2 y3];
        
    elseif(abs(alpha)>0.35)
        k=2/9*(1-alpha);
        h_2=k*(1+2/3*k+7/9*k^2+10/9*k^3+143/81*k^4+728/243*k^5);
        xi_alpha=h_2-2;   %estimate of r1(smallest root)
        % Adusting tje initial estimate: Ostrowksy
        u_0=xi_alpha;
        for iter=1:maxiter
            K = Ostrowksy(u_0,alpha);
            u1=u_0-K;
            u_0=u1;
        end
        r1=u1;
        %phi_alpha=-xi_alpha/2-sqrt(12-3*xi_alpha^2)/2; %estimate of r2 (intermediate root)
        d_plus=(-3*r1+sqrt(-sign(p)*12-3*r1^2))/2;
        d_minus=(-3*r1-sqrt(-sign(p)*12-3*r1^2))/2;
        r2=r1+d_plus;
        r3=r1+d_minus;
        
        %Obtaining the roots of the original cubic
        x1=r1*xm;
        x2=r2*xm;
        x3=r3*xm;
        
        y1=x1-a/3;
        y2=x2-a/3;
        y3=x3-a/3;
        
        y=[y1 y2 y3];
    end  
    
    
    
elseif(p>0)
    %disp('one real roots')
    % Apply the Cardano/Tartaglia formulae
    term1=-1*alpha+sqrt(alpha^2+sign(p));
    term2=-1*alpha-sqrt(alpha^2+sign(p));
    r2=sign(term1)*(abs(term1))^(1/3)+sign(term2)*(abs(term2))^(1/3);
    xm=sqrt(1/3*p);
    x2=r2*xm;
    y2=x2-a/3;
    y=y2;
else
    y=y0;
    %disp('no real roots')
    % Use of fixed-point theory (Newton)
    for iter=1:maxiter
        y=y-(y^3+a*y^2+b*y+c)/(3*y^2+2*a*y+b);
    end
end
end%EOF

function K = Ostrowksy(w,alpha)
K=(w^3-3*w+2*alpha)/(3*(w^2-1))*(1-(2*w*(w^3-3*w+2*alpha))/(3*(w^2-1)))^-1;

end%EOF