function [u]=cmngasa(poly_loc,polyR,Dtg,kpeo,kpmpc,mupeo,mupmpc,i,t,phi_0,theta_0,M,Gpoly,tau0)
format long
%calculate position of source
[x0,y0,z0]=sph2cart(poly_loc(:,5),poly_loc(:,4),polyR);
%poly_loc(:,4)=theta
%poly_loc(:,5)=phi
%convert polymersome position to polar coords
[phip,thetap,rp]=cart2sph(poly_loc(:,1),poly_loc(:,2),poly_loc(:,3));
for n=1:M 
   
e1=@(phi,theta,r,tau) exp(-((r.*sin(phi).*cos(theta)-poly_loc(n,1)).^2+(r.*sin(phi).*sin(theta)-poly_loc(n,2)).^2+(r.*cos(phi)-poly_loc(n,3)).^2)./(4.*Dtg.*(t.*i-tau)));
pterm=@(phi,theta,r,tau) ((4.*pi.*Dtg.*(t.*i-tau)).^(-3/2).*e1(phi,theta,r,tau));

f1=@(phi,theta,tau) mupeo.*(pterm(phi,theta,polyR,tau)-Gpoly).*(kpeo-kpmpc).*(polyR.^2).*sin(phi);
f2=@(phi,theta,r,tau) (mupeo-mupmpc).*(pterm(phi,theta,r,tau)-Gpoly).*(kpmpc).*(r.^2).*sin(phi);
f3=@(phi,theta,r,tau) mupmpc.*(pterm(phi,theta,r,tau)-Gpoly).*kpmpc.*(r.^2).*sin(phi);

I1=integral3(f1,0,phi_0,0,theta_0,0,t*i);
I2=integralN(@(phi,theta,r,tau)f2(phi,theta,r,tau),0,phi_0,0,theta_0,0,polyR,0,t*i);
I3=integralN(@(phi,theta,r,tau)f3(phi,theta,r,tau),0,2*pi,0,pi,0,polyR,0,t*i);
I(n)=I1+I2+I3;
end
%calculate velocity
 u=I'./(4.*pi.*polyR.^2);
 end
