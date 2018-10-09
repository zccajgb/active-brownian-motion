%solution to SS RDE for GDL with area source
%c needs to be integrated over all space for xi eta and zeta


%const mu, non const s
const1=1
const2=2
u=@(phi,theta,phip,thetap,rp) Const1*exp(const2*rp*sin(theta)*sin(thetap)*cos(phi-phip))*exp(const2*cos(theta)*cos(thetap))

integral(@(u,)
