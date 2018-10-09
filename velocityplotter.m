
loops=1
D=0.7e-10;
[k,T,eta,Dtg,petriR,polyR,PDt,PDr,bndclsn,pcolsn,rad,azi,eli,msd]=initialisevariables(M,loops);
[kpmpc,kpeo,theta_0,phi_0]=initialisevariables2;
mupeo=-5e-9;
mupmpc=0;
i=1
timeiter=i;
t=0.005;
tau0=t.*loops;
R=1e-7;

x1=-5e-2:2e-3:5e-2;
 y1=x1;
    
[x,y]=meshgrid(x1,y1);
size(x)
size(y)
pause
poly_loc=[x(:),y(:),zeros(length(x(:)),3)];
M=length(x(:));

deltamu=mupeo-mupmpc;
Gpoly=0;
omega=rotngasa(poly_loc,polyR,Dtg,kpeo,kpmpc,deltamu,timeiter,t,phi_0,theta_0,M,Gpoly,petriR);
u0=cmngasa(poly_loc,polyR,Dtg,kpeo,kpmpc,mupeo,mupmpc,timeiter,t,phi_0,theta_0,M,Gpoly,tau0);
u1=nmngasa(poly_loc,polyR,Dtg,kpeo,kpmpc,deltamu,timeiter,t,phi_0,theta_0,M,Gpoly);
u=u0+u1';
uv=reshape(u,size(x));
omegav=reshape(omega,length(x),length(x));

   a=1
    h=surf(x,-y,-uv)
%      set(h,'edgecolor','interp')
     caxis([0,a/2])
    hold on
     [xa ya za] = cylinder(0.05);
     za(2,:)=a/10;
     j=mesh(xa,ya,za);
     set(j,'edgecolor','black')
    xlabel('x / m')
    ylabel('y / m')
    zlabel('Translational Velocity / m s^{-1}')
    set(gca,'fontsize',10)
    savefig('uvsposn')
    axis([-0.05,0.05,-0.05,0.05,0,a])
clf
h=surf(x,y,omegav)
% axis([x(1),x(end),x(1),x(end),0,1e7])


set(h,'edgecolor','interp')
caxis([0,8e6])
hold on
[xa ya za] = cylinder(0.05);
z(2,:)=1e6;
j=mesh(xa,ya,za);
set(j,'edgecolor','black')
xlabel('x / m')
ylabel('y / m')
zlabel('Difference in G / molec')
set(gca,'fontsize',10)
savefig('wvsposn')