function [omega]=rot_ngpsa(poly_loc,polyR,Dtg,k,mupeo,mupmpc,i,t,phi_0,theta_0,M)
format long
%calculate position of source
[x0,y0,z0]=sph2cart(poly_loc(:,5),poly_loc(:,4),polyR);

%poly_loc(:,4)=theta
%poly_loc(:,5)=phi
%convert polymersome position to polar coords
[phip,thetap,rp]=cart2sph(poly_loc(:,1),poly_loc(:,2),poly_loc(:,3));
for n=1:M 
nhat=[x0(n),y0(n),z0(n)]./norm([x0(n),y0(n),z0(n)]);
n1=nhat(1);
n2=nhat(2);
n3=nhat(3);
%define function to be intergrated, 
% m1a=@(phi,theta,r)(1./r).*(dirac(theta)-dirac(theta-poly_loc(n,4))).*((heaviside(phi)-heaviside(phi-poly_loc(n,5))).*(r.*cos(phi)-z0(n)));
% m1b=@(phi,theta,r)(1./(r.*sin(theta))).*(heaviside(theta)-heaviside(theta-poly_loc(:,4))).*(dirac(phi)-dirac(phi-poly_loc(:,5))).*(r.*sin(theta).*sin(phi)-y0(n));
denom=@(phi,theta,r) ((r.*cos(theta).*sin(phi)-x0(n)).^2+(r.*sin(theta).*sin(phi)-y0(n)).^2+(r.*cos(phi)-z0(n)).^2).^(-3/2);
m1a=@(phi,r) (r.*sin(phi).*(r.*cos(phi)-z0(n)).*n1-(r.*sin(phi).*sin(poly_loc(n,4))-x0(n)).*n3);
m1b=@(phi,r) (r.*sin(phi).*(r.*cos(phi)-z0(n)).*n1-(-x0(n)).*n3);
m2=@(theta,r) (r).*((r.*cos(theta).*sin(poly_loc(n,5))-x0(n)).*n2-(r.*sin(theta).*sin(poly_loc(n,5))-y0(n)).*n1);
f1=@(phi,r) denom(phi,poly_loc(n,4),r).*m1a(phi,r)-denom(phi,0,r).*m1b(phi,r);
f2=@(theta,r) denom(poly_loc(n,5),theta,r).*m2(theta,r);


I1=integral2(f1,0,phi_0,0,polyR);
I2=integral2(f2,0,theta_0,0,polyR);
I(n)=I1+I2;
end
%calculate velocity
e4=exp(-(x0-poly_loc(:,1)).^2+(y0-poly_loc(:,2)).^2+(z0-poly_loc(:,3)).^2);
omega=I'.*3.*(mupeo-mupmpc).*k.*e4./((8.*pi.*polyR^3).*(4.*pi.*Dtg.*t.*i).^(3/2));
end