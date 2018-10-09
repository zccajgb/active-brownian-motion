function [u]=velrotcgasa(x,y,polyR,Dtg,kpeo,kpmpc,deltamu,i,t,phi_0,theta_0,M,Gpoly,petriR)
format long
%calculate position of source
[x0,y0,z0]=sph2cart(0,0,polyR);
%poly_loc(:,4)=theta
%poly_loc(:,5)=phi
%convert polymersome position to polar coords

for n=1:M 
   
%define function to be intergrated, 
% e1=@(phi,theta,r) ((r.*cos(theta).*sin(phi)-x0).^2+(r.*sin(theta).*sin(phi)-y0).^2+(r.*cos(phi)-x0).^2).^(-3/2);
% e2=@(phi,theta,r) (1./r.*(dirac(theta)-dirac(theta-poly_loc(:,4))).*((heaviside(phi)-heaviside(phi-poly_loc(:,5))).*(r.*sin(theta)cos(phi)-y0)
% e3=@(phi,theta,r)(1./r.*(1./(r.^2.*(sin(theta)).^2)).*(heaviside(theta)-heaviside(theta-poly_loc(:,4))).*(dirac(phi)-dirac(phi-poly_loc(:,5))).*(r.*cos(phi)-z0));
nhat=[cos(theta_0).*sin(phi_0);sin(theta_0).*sin(phi_0);cos(phi_0)];
% e1=@(varphi,vartheta,rho) exp(-((rho.*sin(varphi).*cos(vartheta)-poly_loc(n,1)).^2+(rho.*sin(varphi).*sin(vartheta)-poly_loc(n,2)).^2+(rho.*cos(varphi)-poly_loc(n,3)).^2)/(4.*Dtg.*t.*i));
% Denom=@(phi,theta,r,varphi,vartheta,rho) (4.*pi).^(-1).*((r.*cos(theta).*sin(phi)-rho.*cos(vartheta).*sin(varphi)).^2+(r.*sin(theta).*sin(phi)-rho.*sin(vartheta).*sin(varphi)).^2+(r.*cos(phi)-rho.*cos(varphi)).^2).^(-3/2);
% Psi=@(phi,theta,r,varphi,vartheta,rho) ((4.*pi.*Dtg.*t.*i).^(-3/2).*e1(varphi,vartheta,rho)-Gpoly).*Denom(phi,theta,r,varphi,vartheta,rho).*rho.*sin(varphi);
% Psi0=@(phi,theta,r,varphi,vartheta) Psi(phi,theta,r,varphi,vartheta,polyR).*(kpeo-kpmpc);
% Psi1=@(phi,theta,r,varphi,vartheta,rho) Psi(phi,theta,r,varphi,vartheta,rho).*(kpmpc);
% V1=@(phi,theta,r,varphi,vartheta,rho) r.*cos(theta).*sin(phi)-rho.*cos(vartheta).*sin(varphi);
% V2=@(phi,theta,r,varphi,vartheta,rho) r.*sin(theta).*sin(phi)-rho.*sin(vartheta).*sin(varphi);
% V3=@(phi,theta,r,varphi,vartheta,rho) r.*cos(phi)-rho.*cos(varphi);
% V01=@(phi,theta,r,varphi,vartheta) V1(phi,theta,r,varphi,vartheta,polyR);
% V02=@(phi,theta,r,varphi,vartheta) V2(phi,theta,r,varphi,vartheta,polyR);
% V03=@(phi,theta,r,varphi,vartheta) V3(phi,theta,r,varphi,vartheta,polyR);
% V11=@(phi,theta,r,varphi,vartheta,rho) V1(phi,theta,r,varphi,vartheta,rho);
% V12=@(phi,theta,r,varphi,vartheta,rho) V2(phi,theta,r,varphi,vartheta,rho);
% V13=@(phi,theta,r,varphi,vartheta,rho) V3(phi,theta,r,varphi,vartheta,rho);
% 
% f1=@(phi,r,varphi,vartheta) Psi0(phi,0,r,varphi,vartheta).*V02(phi,0,r,varphi,vartheta).*r.*sin(phi).*nhat(1);
% f2=@(phi,r,varphi,vartheta) Psi0(phi,theta_0,r,varphi,vartheta).*V02(phi,theta_0,r,varphi,vartheta).*r.*sin(phi);
% f3=@(theta,r,varphi,vartheta) Psi0(0,theta,r,varphi,vartheta).*V03(0,theta,r,varphi,vartheta).*r.*nhat(1);
% f4=@(theta,r,varphi,vartheta) -Psi0(phi_0,theta,r,varphi,vartheta).*V03(phi_0,theta,r,varphi,vartheta).*r.*nhat(1);
% f5=@(theta,r,varphi,vartheta) -Psi0(0,theta,r,varphi,vartheta).*V01(0,theta,r,varphi,vartheta).*r.*nhat(2);
% f6=@(theta,r,varphi,vartheta) -Psi0(phi_0,theta,r,varphi,vartheta).*V01(phi_0,theta,r,varphi,vartheta).*r.*nhat(2);
% f7=@(phi,r,varphi,vartheta) -Psi0(phi,0,r,varphi,vartheta).*V01(phi,0,r,varphi,vartheta).*r.*sin(phi).*nhat(3);
% f8=@(phi,r,varphi,vartheta) Psi0(phi,theta_0,r,varphi,vartheta).*V01(phi,theta_0,r,varphi,vartheta).*r.*sin(phi).*nhat(3);
% 
% f9=@(phi,r,varphi,vartheta,rho) Psi1(phi,0,r,varphi,vartheta,rho).*V12(phi,0,r,varphi,vartheta,rho).*r.*sin(phi).*nhat(1);
% f10=@(phi,r,varphi,vartheta,rho) Psi1(phi,theta_0,r,varphi,vartheta,rho).*V12(phi,theta_0,r,varphi,vartheta,rho).*r.*sin(phi);
% f11=@(theta,r,varphi,vartheta,rho) Psi1(0,theta,r,varphi,vartheta,rho).*V13(0,theta,r,varphi,vartheta,rho).*r.*nhat(1);
% f12=@(theta,r,varphi,vartheta,rho) -Psi1(phi_0,theta,r,varphi,vartheta,rho).*V13(phi_0,theta,r,varphi,vartheta,rho).*r.*nhat(1);
% f13=@(theta,r,varphi,vartheta,rho) -Psi1(0,theta,r,varphi,vartheta,rho).*V11(0,theta,r,varphi,vartheta,rho).*r.*nhat(2);
% f14=@(theta,r,varphi,vartheta,rho) -Psi1(phi_0,theta,r,varphi,vartheta,rho).*V11(phi_0,theta,r,varphi,vartheta,rho).*r.*nhat(2);
% f15=@(phi,r,varphi,vartheta,rho) -Psi1(phi,0,r,varphi,vartheta,rho).*V11(phi,0,r,varphi,vartheta,rho).*r.*sin(phi).*nhat(3);
% f16=@(phi,r,varphi,vartheta,rho) Psi1(phi,theta_0,r,varphi,vartheta,rho).*V11(phi,theta_0,r,varphi,vartheta,rho).*r.*sin(phi).*nhat(3);
% 
% F1=@(phi,r,varphi,vartheta) (kpeo-kpmpc).*(f1(phi,r,varphi,vartheta)+f2(phi,r,varphi,vartheta)+f7(phi,r,varphi,vartheta)+f8(phi,r,varphi,vartheta)); 
% F2=@(theta,r,varphi,vartheta) (kpeo-kpmpc).*(f3(theta,r,varphi,vartheta)+f4(theta,r,varphi,vartheta)+f5(theta,r,varphi,vartheta)+f6(theta,r,varphi,vartheta));
% F3=@(phi,r,varphi,vartheta,rho) kpmpc.*(f9(phi,r,varphi,vartheta,rho)+f10(phi,r,varphi,vartheta,rho)+f15(phi,r,varphi,vartheta,rho)+f16(phi,r,varphi,vartheta,rho));
% F4=@(theta,r,varphi,vartheta,rho) kpmpc.*(f11(theta,r,varphi,vartheta,rho)+f12(theta,r,varphi,vartheta,rho)+f13(theta,r,varphi,vartheta,rho)+f14(theta,r,varphi,vartheta,rho));

I1=montecarlointr(@rF1,[0,0,0,0],[phi_0,polyR,phi_0,theta_0],1e-3,1,0,poly_loc(n,:),Dtg,t,i,Gpoly,kpeo,kpmpc,polyR,phi_0,theta_0,nhat);
I2=montecarlointr(@rF2,[0,0,0,0],[theta_0,polyR,phi_0,theta_0],1e-3,0,1,poly_loc(n,:),Dtg,t,i,Gpoly,kpeo,kpmpc,polyR,phi_0,theta_0,nhat);
I3=montecarlointr(@rF3,[0,0,0,0,0],[phi_0,polyR,phi_0,theta_0,polyR],1e-3,1,0,poly_loc(n,:),Dtg,t,i,Gpoly,kpeo,kpmpc,polyR,phi_0,theta_0,nhat);
I4=montecarlointr(@rF4,[0,0,0,0,0],[theta_0,polyR,phi_0,theta_0,polyR],1e-3,0,1,poly_loc(n,:),Dtg,t,i,Gpoly,kpeo,kpmpc,polyR,phi_0,theta_0,nhat);

% I4=integralN(@(theta,r,varphi,vartheta,rho)F4(theta,r,varphi,vartheta,rho),0,theta_0,0,polyR,0,phi_0,0,theta_0,0,petriR,'method','iterated');
I(n)=I1+I2+I3+I4;

end
%calculate velocity
u=3.*deltamu.*I./(2.*4.*pi.*polyR.^3);
end
