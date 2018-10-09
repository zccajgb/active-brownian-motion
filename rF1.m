function [F1]=rF1(phi,theta_0,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly,kpeo,kpmpc,nhat)
f1=uPsi(phi,0,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly).*uV2(phi,0,r,varphi,vartheta,polyR).*r.*sin(phi).*nhat(1);
f2=uPsi(phi,theta_0,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly).*uV2(phi,theta_0,r,varphi,vartheta,polyR).*r.*sin(phi).*nhat(1);
f7=-uPsi(phi,0,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly).*uV1(phi,0,r,varphi,vartheta,polyR).*r.*sin(phi).*nhat(3);
f8=uPsi(phi,theta_0,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly).*uV1(phi,theta_0,r,varphi,vartheta,polyR).*r.*sin(phi).*nhat(3);

F1=(kpeo-kpmpc).*(f1+f2+f7+f8);
end