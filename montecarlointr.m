function [integrand]=montecarlointr(f,l,u,abstol,phival,thetaval,pln,Dtg,t,i,Gpoly,kpeo,kpmpc,polyR,phi0,theta0,nhat)
N=1e5;
n=N;
m=length(u);
rnd=rand(N,m);
x=(u-l).*rnd+l;
V=prod(u-l);
if phival==0
    if (m==4)
    I=f(phi0,x(:,1),x(:,2),x(:,3),x(:,4),polyR,pln,Dtg,t,i,Gpoly,kpeo,kpmpc,nhat);
    elseif (m==5)
    I=f(phi0,x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),pln,Dtg,t,i,Gpoly,kpeo,kpmpc,nhat);   
    end
elseif thetaval==0
    if (m==4)
    I=f(x(:,1),theta0,x(:,2),x(:,3),x(:,4),polyR,pln,Dtg,t,i,Gpoly,kpeo,kpmpc,nhat);
    elseif (m==5)
    I=f(x(:,1),theta0,x(:,2),x(:,3),x(:,4),x(:,5),pln,Dtg,t,i,Gpoly,kpeo,kpmpc,nhat);   
    end  
end
    mn=mean(I);
    integrand=V.*mn;
    abserr=N.^(-1/2);
   index=1;
   
% while abs(abserr)>abstol
%   rnd=rand(n,length(u));
%   x=(u-l).*rnd+l;
%   N=N+n;
%   if m==4
%     I2=f(x(:,1),x(:,2),x(:,3),x(:,4));
%   end    
%   if m==5
%     I2=f(x(:,1),x(:,2),x(:,3),x(:,4),x(:,5));
%   end
%   I=[I;I2];
%     mn=mean(I);
%     integrand=V.*mn;
%    abserr=N.^(-1/2);
% index=index+1;
% end
end
