function [kpmpc,kpeo,theta_0,phi_0]=initialisevariablesGs2
Dpmpc=1e-3;
dpmpc=6.4e-9;
dpeo=2.4e-9;
kpmpc=Dpmpc/dpmpc;
Dpeo=0.1;
kpeo=Dpeo/dpeo;
Rpeo=15e-9;
Rpmpc=50e-9;
theta_0=2.*asin(Rpeo./Rpmpc);
phi_0=theta_0;
end