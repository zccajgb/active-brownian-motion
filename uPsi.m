function [Psi]=uPsi(phi,theta,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly)
xi=rho.*sin(varphi).*cos(vartheta);
eta=rho.*sin(varphi).*sin(vartheta);
zeta=rho.*sin(varphi).*sin(vartheta);
sp=sin(phi);
% e1=exp(-((xi-pln(1)).^2+(eta-pln(2)).^2+(zeta-pln(3)).^2)/(4.*Dtg.*10000));
 powa=(r.*cos(theta).*sp-xi).^2+(r.*sin(theta).*sp-eta).^2+(r.*cos(phi)-zeta).^2;
rsq=realsqrt(powa);
Denom=(4.*pi).*rsq.*rsq.*rsq;
 G=(4.*pi).*rsq;
% Psi=((4.*pi.*Dtg).^(-3/2).*e1-Gpoly).*rho.*sin(varphi)./Denom;
 Psi=((1./G)-Gpoly).*rho.*sin(varphi)./Denom;
end
