function [V2]=uV2(phi,theta,r,varphi,vartheta,rho)
V2=r.*sin(theta).*sin(phi)-rho.*sin(vartheta).*sin(varphi);
end