function [F3]=rF3(phi,theta_0,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly,kpeo,kpmpc,nhat)
f9=uPsi(phi,0,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly).*uV2(phi,0,r,varphi,vartheta,rho).*r.*sin(phi).*nhat(1);
f10=uPsi(phi,theta_0,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly).*uV2(phi,theta_0,r,varphi,vartheta,rho).*r.*sin(phi);
f15=-uPsi(phi,0,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly).*uV1(phi,0,r,varphi,vartheta,rho).*r.*sin(phi).*nhat(3);
f16=uPsi(phi,theta_0,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly).*uV1(phi,theta_0,r,varphi,vartheta,rho).*r.*sin(phi).*nhat(3);
F3=kpmpc.*(f9+f10+f15+f16);
end