function [ x,flag,fun_p_val ] = Cubic_rootsv2( a_tilde,b_tilde,c_tilde,d_tilde,mu)
%Cubic_roots returns the real positive roots of 
%the polynomial of order 3 in the normal form:
%x^3+px^2+qx+r
%%%%Inputs%%%%
% a_tilde: coefficient of the cubic term
% b_tilde: coefficient of the quadratic term
% c_tilde: coefficient of the first-order term
% d_tilde: constant term
% mu:Lagrangian Multiplier
%%%%Outputs%%%
% x: roots
% fun_p: evaluation of the derivative of the root at mu
% flag : 1 if only one real root, o otherwise
%% Computation of the canonical form:x^3+px^2+qx+r
% p: coefficient of the quadratic term
% q: coefficient of the first-order term
% r: constant term
p=(b_tilde+mu)/a_tilde;q=c_tilde/a_tilde;r=(d_tilde)/a_tilde;

%% Computation of the normal form y^3+ay+b=0 with following CV x=y-p/3
a=1/3*(3*q-p^2);
b=1/27*(2*p^3-9*p*q+27*r);
flag=0;
%% Computation and discussion based on the radicants
rad=b^2/4+a^3/27;
y=zeros(3,1);
if(rad>0)
    %fprintf('One real root\n');
    A=nthroot(-b/2+sqrt(rad),3);
    B=nthroot(-b/2-sqrt(rad),3);
    y(1)=A+B;
    y(2)=-1/2*(A+B)+1i*sqrt(3)/2*(A-B);
    y(3)=-1/2*(A+B)-1i*sqrt(3)/2*(A-B);
    %Computation of the derivitative of f(mu_k)=sum_i w_i-c
    fac1=1/(nthroot((-b/2+sqrt(rad))^2,3)+eps);
    fac2=1/2*1/(sqrt(rad)+eps)*d_tilde/(a_tilde+eps)-1;
    fac3=1/(nthroot((-b/2-sqrt(rad))^2,3)+eps);
    fac4=-1/2*1/(sqrt(rad)+eps)*d_tilde/(a_tilde+eps)-1;
    fun_p_val=1/(27*a_tilde^3+eps)*(b_tilde+mu)^2*(fac1*fac2+fac3*fac4)-1/3*1/(a_tilde+eps);
    flag=1;
elseif(rad==0)
    %fprintf('Three real roots of which at least two are equal\n');
    if(b>0)
        y(1)=-2*sqrt(-a/3);
        y(2)=sqrt(-a/3);
        y(3)=sqrt(-a/3);
    elseif(b<0)
        y(1)=2*sqrt(-a/3);
        y(2)=-1*sqrt(-a/3);
        y(3)=-1*sqrt(-a/3);
    else
        y=zeros(3,1);
    end
    
else
%     syms muk;
%     psym=(b_tilde+muk)/a_tilde;
%     asym=1/3*(3*q-psym^2);
%     bsym=1/27*(2*psym^3-9*p*q+27*r);
    %fprintf('Three real and unequal roots\n');
    s=b_tilde+mu;
    t=(d_tilde/(a_tilde)+2/27*(s)^3/(a_tilde^3));
    if(b>0)
        phi=acos(-1*sqrt((b^2/4)/(-a^3/27)));
        %phi_sym=acos(-1*sqrt((bsym^2/4)/(-asym^3/27)));
        num=2*(b_tilde+mu)*sin(1/3*acos(27/2*sqrt((a_tilde^6*t^2)/(s)^6))+pi/6);
        denom=3*a_tilde^2*sqrt((s^2)/(a_tilde^2));
        A=(num)/(denom+eps);
        B=(3*sqrt((s^2)/(a_tilde^2+eps))*(((4*a_tilde^3*t)/(9*s^4+eps))-((6*a_tilde^6*t^2)/(s^7+eps))))*cos(1/3*acos(27/2*sqrt((a_tilde^6*t^2)/(s)^6))+pi/6);
        C=2*sqrt((a_tilde^6*t^2)/(s^6+eps))*sqrt(1-(729*a_tilde^6*t^2)/(4*s^6+eps));
        fun_p_val(1)=A-B/C-1/3*1/(a_tilde+eps);
    elseif(b<0)
        phi=acos(sqrt((b^2/4)/(-a^3/27)));
        %phi_sym=acos(sqrt((bsym^2/4)/(-asym^3/27)));
        num=2*(b_tilde+mu)*cos(1/3*acos(27/2*sqrt((a_tilde^6*t^2)/(s)^6)));
        denom=3*a_tilde^2*sqrt((s^2)/(a_tilde^2+eps));
        A=(num)/(denom+eps);
        B=(3*sqrt((s^2)/(a_tilde^2+eps))*(((4*a_tilde^3*t)/(9*s^4+eps))-((6*a_tilde^6*t^2)/(s^7+eps))))*sin(1/3*acos(27/2*sqrt((a_tilde^6*t^2)/(s)^6)));
        C=2*sqrt((a_tilde^6*t^2)/(s^6+eps))*sqrt(1-(729*a_tilde^6*t^2)/(4*s^6+eps));
        fun_p_val(1)=A+B/C-1/3*1/(a_tilde+eps);
    end
    k=0;
    y(1)=2*sqrt(-a/3)*cos(phi/3+2*k*pi/3);
%     for i=1:3
%         k=i-1;
%         y(i)=2*sqrt(-a/3)*cos(phi/3+2*k*pi/3);
%         ysym(i)=2*sqrt(-asym/3)*cos(phi_sym/3+2*k*pi/3)-psym/3;
%         fun_p(i)=diff(ysym(i),muk);
%         fun_p_val(i)=eval(subs(fun_p(i),[muk],[mu]));
%     end
    
end
x=y-p/3;
end
