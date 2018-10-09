function [V3]=uV3(phi,theta,r,varphi,vartheta,rho)
V3=r.*cos(phi)-rho.*cos(varphi);
end