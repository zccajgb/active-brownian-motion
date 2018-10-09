function [F4]=uF4(phi_0,theta,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly,kpeo,kpmpc)
f7=uPsi(0,theta,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly).*(kpmpc).*uV3(0,theta,r,varphi,vartheta,rho).*r;
f8=-uPsi(phi_0,theta,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly).*(kpmpc).*uV3(phi_0,theta,r,varphi,vartheta,rho).*r;
F4=f7+f8;
end