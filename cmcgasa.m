function [u]=cmcgasa(poly_loc,polyR,Dtg,kpeo,kpmpc,mupeo,mupmpc,i,t,phi_0,theta_0,M,Gpoly,tau0)
format long
%calculate position of source
[x0,y0,z0]=sph2cart(poly_loc(:,5),poly_loc(:,4),polyR);
%poly_loc(:,4)=theta
%poly_loc(:,5)=phi
%convert polymersome position to polar coords
[phip,thetap,rp]=cart2sph(poly_loc(:,1),poly_loc(:,2),poly_loc(:,3));
for n=1:M 
  
%  e1=@(phi,theta,r) exp(-((r.*sin(phi).*cos(theta)-poly_loc(n,1)).^2+(r.*sin(phi).*sin(theta)-poly_loc(n,2)).^2+(r.*cos(phi)-poly_loc(n,3)).^2)./(4.*10000.*Dtg));
%  pterm=@(phi,theta,r) ((4.*pi.*Dtg).^(-3/2).*e1(phi,theta,r));
pterm=@(phi,theta,r) 1./sqrt((r.*sin(phi).*cos(theta)-poly_loc(n,1)).^2+(r.*sin(phi).*sin(theta)-poly_loc(n,2)).^2+(r.*cos(phi)-poly_loc(n,3)).^2);
 
f1=@(phi,theta) mupeo.*(pterm(phi,theta,polyR)-Gpoly).*(kpeo-kpmpc).*(polyR.^2).*sin(phi);
f2=@(phi,theta,r) (mupeo-mupmpc).*(pterm(phi,theta,r)-Gpoly).*(kpmpc).*(r.^2).*sin(phi);
f3=@(phi,theta,r) mupmpc.*(pterm(phi,theta,r)-Gpoly).*kpmpc.*(r.^2).*sin(phi);

I1=integral2(f1,0,phi_0,0,theta_0);
I2=integral3(@(phi,theta,r)f2(phi,theta,r),0,phi_0,0,theta_0,0,polyR);
I3=integral3(@(phi,theta,r)f3(phi,theta,r),0,2*pi,0,pi,0,polyR);
I(n)=I1+I2+I3;
end
%calculate velocity
 u=I'./(4.*pi.*polyR.^2);
 end
