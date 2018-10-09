function [F4]=rF4(phi_0,theta,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly,kpeo,kpmpc,nhat)
f11=uPsi(0,theta,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly).*uV3(0,theta,r,varphi,vartheta,rho).*r.*nhat(1);
f12=-uPsi(phi_0,theta,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly).*uV3(phi_0,theta,r,varphi,vartheta,rho).*r.*nhat(1);
f13=-uPsi(0,theta,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly).*uV1(0,theta,r,varphi,vartheta,rho).*r.*nhat(2);
f14=-uPsi(phi_0,theta,r,varphi,vartheta,rho,pln,Dtg,t,i,Gpoly).*uV1(phi_0,theta,r,varphi,vartheta,rho).*r.*nhat(2);
F4=kpmpc.*(f11+f12+f13+f14);
end