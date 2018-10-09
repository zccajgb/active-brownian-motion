function [out]=exponfunc(r,phi,theta,poly_loc,Dtg,t,i,n)
out=exp(-((r.*sin(phi).*cos(theta)-poly_loc(n,1)).^2+(r.*sin(phi).*sin(theta)-poly_loc(n,2)).^2+(r.*cos(phi)-poly_loc(n,3)).^2)/(4.*Dtg.*t.*i))
end