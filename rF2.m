function [F2]=rF2(phi_0,theta,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly,kpeo,kpmpc,nhat)
f3=uPsi(0,theta,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly).*uV3(0,theta,r,varphi,vartheta,polyR).*r.*nhat(1);
f4=-uPsi(phi_0,theta,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly).*uV3(phi_0,theta,r,varphi,vartheta,polyR).*r.*nhat(1);
f5=-uPsi(0,theta,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly).*uV1(0,theta,r,varphi,vartheta,polyR).*r.*nhat(2);
f6=-uPsi(phi_0,theta,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly).*uV1(phi_0,theta,r,varphi,vartheta,polyR).*r.*nhat(2);
F2=(kpeo-kpmpc).*(f3+f4+f5+f6);
end