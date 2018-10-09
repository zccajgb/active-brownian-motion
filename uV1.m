function [V1]=uV1(phi,theta,r,varphi,vartheta,rho)
V1=r.*cos(theta).*sin(phi)-rho.*cos(vartheta).*sin(varphi);
end