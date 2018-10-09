function [u]=nmngpsa(poly_loc,polyR,Dtg,k,mupeo,mupmpc,i,t,phi_0,theta_0,M)
format long
%calculate position of source
[x0,y0,z0]=sph2cart(poly_loc(:,5),poly_loc(:,4),polyR);
%poly_loc(:,4)=theta
%poly_loc(:,5)=phi
%convert polymersome position to polar coords
[phip,thetap,rp]=cart2sph(poly_loc(:,1),poly_loc(:,2),poly_loc(:,3));
for n=1:M 
   
%define function to be intergrated, 
% e1=@(phi,theta,r) ((r.*cos(theta).*sin(phi)-x0).^2+(r.*sin(theta).*sin(phi)-y0).^2+(r.*cos(phi)-x0).^2).^(-3/2);
% e2=@(phi,theta,r) (1./r.*(dirac(theta)-dirac(theta-poly_loc(:,4))).*((heaviside(phi)-heaviside(phi-poly_loc(:,5))).*(rsin(theta)cos(phi)-y0)
% e3=@(phi,theta,r)(1./r.*(1./(r.^2.*(sin(theta)).^2)).*(heaviside(theta)-heaviside(theta-poly_loc(:,4))).*(dirac(phi)-dirac(phi-poly_loc(:,5))).*(r.*cos(phi)-z0));

% f=@(phi,theta) e1.*(e2+e3));
denom=@(phi,theta,r) ((r.*cos(theta).*sin(phi)-x0(n)).^2+(r.*sin(theta).*sin(phi)-y0(n)).^2+(r.*cos(phi)-z0(n)).^2).^(-3/2);
m1a=@(phi,r) (r.*sin(poly_loc(n,4)).*sin(phi)-y0(n)).*r.*sin(phi);
m1b=@(phi,r) (y0(n)).*r.*sin(phi);

m2=@(theta,r) (r.*cos(poly_loc(n,5))-z0(n)).*r;
f1=@(phi,r) denom(phi,poly_loc(n,4),r).*m1a(phi,r)+denom(phi,0,r).*m1b(phi,r);
f2=@(theta,r) denom(poly_loc(n,5),theta,r).*m2(theta,r);
I1=integral2(f1,0,phi_0,0,polyR);
I2=integral2(f1,0,theta_0,0,polyR);
I(n)=I1+I2;

end
%calculate velocity
e4=exp(-(x0-poly_loc(:,1)).^2+(y0-poly_loc(:,2)).^2+(z0-poly_loc(:,3)).^2);
u=k.*(mupeo-mupmpc).*e4./(4.*pi.*(4.*pi.*Dtg.*t.*i).^(3/2)).*I';
end

