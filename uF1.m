function [F1]=uF1(phi,theta_0,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly,kpeo,kpmpc)
f1=uPsi(phi,0,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly).*(kpeo-kpmpc).*uV2(phi,0,r,varphi,vartheta,polyR).*r.*sin(phi);
f2=-uPsi(phi,theta_0,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly).*(kpeo-kpmpc).*uV2(phi,theta_0,r,varphi,vartheta,polyR).*r.*sin(phi);
F1=f1+f2;
end