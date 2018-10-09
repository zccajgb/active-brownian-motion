function [u]=nmngasa(poly_loc,polyR,Dtg,kpeo,kpmpc,deltamu,i,t,phi_0,theta_0,M,Gpoly)
format long
%calculate position of source
[x0,y0,z0]=sph2cart(poly_loc(:,5),poly_loc(:,4),polyR);
%poly_loc(:,4)=theta
%poly_loc(:,5)=phi
%convert polymersome position to polar coords
[phip,thetap,rp]=cart2sph(poly_loc(:,1),poly_loc(:,2),poly_loc(:,3));
for n=1:M 
   
%define function to be intergrated, 
e1=@(phi,theta,r) ((r.*cos(theta).*sin(phi)-x0).^2+(r.*sin(theta).*sin(phi)-y0).^2+(r.*cos(phi)-x0).^2).^(-3/2);
e2=@(phi,theta,r) (1./r.*(dirac(theta)-dirac(theta-poly_loc(:,4))).*((heaviside(phi)-heaviside(phi-poly_loc(:,5))).*(r.*sin(theta).*cos(phi)-y0)));
e3=@(phi,theta,r)(1./r.*(1./(r.^2.*(sin(theta)).^2)).*(heaviside(theta)-heaviside(theta-poly_loc(:,4))).*(dirac(phi)-dirac(phi-poly_loc(:,5))).*(r.*cos(phi)-z0));
e1=@(varphi,vartheta,rho) exp(-((rho.*sin(varphi).*cos(vartheta)-poly_loc(n,1)).^2+(rho.*sin(varphi).*sin(vartheta)-poly_loc(n,2)).^2+(rho.*cos(varphi)-poly_loc(n,3)).^2)/(4.*Dtg.*t.*i));
Denom=@(phi,theta,r,varphi,vartheta,rho) (4.*pi).^(-1).*((r.*cos(theta).*sin(phi)-rho.*cos(vartheta).*sin(varphi)).^2+(r.*sin(theta).*sin(phi)-rho.*sin(vartheta).*sin(varphi)).^2+(r.*cos(phi)-rho.*cos(varphi)).^2).^(-3/2);
Psi=@(phi,theta,r,varphi,vartheta,rho) ((4.*pi.*Dtg.*t.*i).^(-3/2).*e1(varphi,vartheta,rho)-Gpoly).*Denom(phi,theta,r,varphi,vartheta,rho).*rho.*sin(varphi);
Psi0=@(phi,theta,r,varphi,vartheta) Psi(phi,theta,r,varphi,vartheta,polyR).*(kpeo-kpmpc);
Psi1=@(phi,theta,r,varphi,vartheta,rho) Psi(phi,theta,r,varphi,vartheta,rho).*(kpmpc);
V1=@(phi,theta,r,varphi,vartheta,rho) r.*cos(theta).*sin(phi)-rho.*cos(vartheta).*sin(varphi);
V2=@(phi,theta,r,varphi,vartheta,rho) r.*sin(theta).*sin(phi)-rho.*sin(vartheta).*sin(varphi);
V3=@(phi,theta,r,varphi,vartheta,rho) r.*cos(phi)-rho.*cos(varphi);
V01=@(phi,theta,r,varphi,vartheta) V1(phi,theta,r,varphi,vartheta,polyR);
V02=@(phi,theta,r,varphi,vartheta) V2(phi,theta,r,varphi,vartheta,polyR);
V03=@(phi,theta,r,varphi,vartheta) V3(phi,theta,r,varphi,vartheta,polyR);
V11=@(phi,theta,r,varphi,vartheta,rho) V1(phi,theta,r,varphi,vartheta,rho);
V12=@(phi,theta,r,varphi,vartheta,rho) V2(phi,theta,r,varphi,vartheta,rho);
V13=@(phi,theta,r,varphi,vartheta,rho) V3(phi,theta,r,varphi,vartheta,rho);

f1=@(phi,r,varphi,vartheta,polyR) Psi0(phi,theta_0,r,varphi,vartheta).*V02(phi,0,r,varphi,vartheta).*r.*sin(phi);
f2=@(phi,r,varphi,vartheta,polyR) -Psi0(phi,theta_0,r,varphi,vartheta).*V02(phi,theta_0,r,varphi,vartheta).*r.*sin(phi);
f3=@(theta,r,varphi,vartheta,polyR) Psi0(phi_0,theta,r,varphi,vartheta).*V03(0,theta,r,varphi,vartheta).*r;
f4=@(theta,r,varphi,vartheta,polyR) -Psi0(phi_0,theta,r,varphi,vartheta).*V03(phi_0,theta,r,varphi,vartheta).*r;
f5=@(phi,r,varphi,vartheta,rho) Psi1(phi,theta_0,r,varphi,vartheta,rho).*V12(phi,0,r,varphi,vartheta,rho).*r.*sin(phi);
f6=@(phi,r,varphi,vartheta,rho) -Psi1(phi,theta_0,r,varphi,vartheta,rho).*(kpmpc).*V12(phi,theta_0,r,varphi,vartheta,rho).*r.*sin(phi);
f7=@(theta,r,varphi,vartheta,rho) Psi1(phi_0,theta,r,varphi,vartheta,rho).*(kpmpc).*V13(0,theta,r,varphi,vartheta,rho).*r;
f8=@(theta,r,varphi,vartheta,rho) -Psi1(phi_0,theta,r,varphi,vartheta,rho).*(kpmpc).*V13(phi_0,theta,r,varphi,vartheta,rho).*r;

F1=@(phi,r,varphi,vartheta) (f1(phi,r,varphi,vartheta)+f2(phi,r,varphi,vartheta ));
F2=@(theta,r,varphi,vartheta) (f3(theta,r,varphi,vartheta)+f4(theta,r,varphi,vartheta));
F3=@(phi,r,varphi,vartheta,rho) (f5(phi,r,varphi,vartheta,rho)+f6(phi,r,varphi,vartheta,rho));
F4=@(theta,r,varphi,vartheta,rho) (f7(theta,r,varphi,vartheta,rho)+f8(theta,r,varphi,vartheta,rho));
% tic
% I1=montecarloint(@uF1,[0,0,0,0],[phi_0,polyR,phi_0,theta_0],1e-3,1,0,poly_loc(n,:),Dtg,t,i,Gpoly,kpeo,kpmpc,polyR,phi_0,theta_0);
% I2=montecarloint(@uF2,[0,0,0,0],[theta_0,polyR,phi_0,theta_0],1e-3,0,1,poly_loc(n,:),Dtg,t,i,Gpoly,kpeo,kpmpc,polyR,phi_0,theta_0);
% I3=montecarloint(@uF3,[0,0,0,0,0],[phi_0,polyR,phi_0,theta_0,polyR],1e-3,1,0,poly_loc(n,:),Dtg,t,i,Gpoly,kpeo,kpmpc,polyR,phi_0,theta_0);
% I4=montecarloint(@uF4,[0,0,0,0,0],[theta_0,polyR,phi_0,theta_0,polyR],1e-3,0,1,poly_loc(n,:),Dtg,t,i,Gpoly,kpeo,kpmpc,polyR,phi_0,theta_0);
% toc
I1=integralN(@(phi,r,varphi,vartheta)F1(phi,r,varphi,vartheta),0,phi_0,0,polyR,0,phi_0,0,theta_0);
I2=integralN(@(theta,r,varphi,vartheta)F2(theta,r,varphi,vartheta),0,theta_0,0,polyR,0,phi_0,0,theta_0);
I3=integralN(@(phi,r,varphi,vartheta,rho)F3(phi,r,varphi,vartheta,rho),0,phi_0,0,polyR,0,phi_0,0,theta_0,0,polyR);
I4=integralN(@(theta,r,varphi,vartheta,rho)F4(theta,r,varphi,vartheta,rho),0,theta_0,0,polyR,0,phi_0,0,theta_0,0,polyR);

I(n)=I1+I2+I3+I4;

end
%calculate velocity
u=deltamu.*I./polyR;
end
