function [u]=velcmcgasa(x,y,polyR,Dtg,kpeo,kpmpc,mupeo,mupmpc,i,t,phi_0,theta_0,M,Gpoly,tau0)
format long
%calculate position of source
[x0,y0]=pol2cart(0,polyR);
%poly_loc(:,4)=theta
%poly_loc(:,5)=phi
%convert polymersome position to polar coords

for n=1:M 
  
%  e1=@(phi,theta,r) exp(-((r.*sin(phi).*cos(theta)-poly_loc(n,1)).^2+(r.*sin(phi).*sin(theta)-poly_loc(n,2)).^2+(r.*cos(phi)-poly_loc(n,3)).^2)./(4.*10000.*Dtg));
%  pterm=@(phi,theta,r) ((4.*pi.*Dtg).^(-3/2).*e1(phi,theta,r));
pterm=@(theta,r) 1./sqrt((r*cos(theta)-poly_loc(n,1)).^2+(r.*sin(theta)-poly_loc(n,2)).^2);
 
f1=@(theta) mupeo.*(pterm(theta,polyR)-Gpoly).*(kpeo-kpmpc).*(polyR.^2);
f2=@(theta,r) (mupeo-mupmpc).*(pterm(theta,r)-Gpoly).*(kpmpc).*(r.^2);
f3=@(theta,r) mupmpc.*(pterm(theta,r)-Gpoly).*kpmpc.*(r.^2);

I1=integral1(f1,0,theta_0);
I2=integral2(@(theta,r)f2(theta,r),0,theta_0,0,polyR);
I3=integral3(@(phi,theta,r)f3(phi,theta,r),0,2*pi,0,pi,0,polyR);
I(n)=I1+I2+I3;
end
%calculate velocity
 u=I'./(4.*pi.*polyR.^2);
 end
