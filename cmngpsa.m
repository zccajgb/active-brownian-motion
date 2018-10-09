function [u]=cmngpsa(poly_loc,polyR,Dtg,k,mu,i,t)
[x0,y0,z0]=sph2cart(poly_loc(:,5),poly_loc(:,4),polyR);
ro=((poly_loc(:,1)-x0).^2+(poly_loc(:,2)-x0).^2+(poly_loc(:,3)-z0).^2);
u=(k*mu/(8*pi*Dtg))*(4*pi*Dtg*i*t)^(3/2).*exp(-ro./(4*Dtg*t*i));
end
