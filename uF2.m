function [F2]=uF2(phi_0,theta,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly,kpeo,kpmpc)
f3=uPsi(0,theta,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly).*(kpeo-kpmpc).*uV3(0,theta,r,varphi,vartheta,polyR).*r;
f4=-uPsi(phi_0,theta,r,varphi,vartheta,polyR,pln,Dtg,t,i,Gpoly).*(kpeo-kpmpc).*uV3(phi_0,theta,r,varphi,vartheta,polyR).*r;
F2=f3+f4;
end