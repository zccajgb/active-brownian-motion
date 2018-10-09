function [F3]=uF3(phi,theta_0,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly,kpeo,kpmpc)
f5=uPsi(phi,0,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly).*(kpmpc).*uV2(phi,0,r,varphi,vartheta,rho).*r.*sin(phi);
f6=-uPsi(phi,theta_0,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly).*(kpmpc).*uV2(phi,theta_0,r,varphi,vartheta,rho).*r.*sin(phi);
F3=f5+f6;
end