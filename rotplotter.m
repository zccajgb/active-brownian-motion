clear all
xl=[-0.04,-0.0001];
for ind=1:2
loops=1
D=0.7e-10;
theta=[0.1:0.1:2*pi]';
M=length(theta);
[k,T,eta,Dtg,petriR,polyR,PDt,PDr,bndclsn,pcolsn,rad,azi,eli,msd]=initialisevariables(M,loops);
[kpmpc,kpeo,theta_0,phi_0]=initialisevariables2;
mupeo=-5e-9;
mupmpc=0;
i=200;
timeiter=i;
t=0.005;
tau0=t.*loops;
x=xl(ind);
poly_loc=x[x*ones(length(theta),1),zeros(length(theta),2),theta,zeros(length(theta),1)];


deltamu=mupeo-mupmpc;
Gpoly=0;
omega=rotngasa(poly_loc,polyR,Dtg,kpeo,kpmpc,deltamu,timeiter,t,phi_0,theta_0,M,Gpoly,petriR);
u0=cmngasa(poly_loc,polyR,Dtg,kpeo,kpmpc,mupeo,mupmpc,timeiter,t,phi_0,theta_0,M,Gpoly,tau0);
u1=nmngasa(poly_loc,polyR,Dtg,kpeo,kpmpc,deltamu,timeiter,t,phi_0,theta_0,M,Gpoly);
u=u0+u1';
uv(:,ind)=u;
omegav(:,ind)=omega;
end