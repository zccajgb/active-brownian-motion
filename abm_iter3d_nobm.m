function [poly_loc]= abm_iter3d_nobm(PDr,PDt,t,n,poly_loc,u,omega); 

U=[u.*sin(poly_loc(:,5)).*cos(poly_loc(:,4)),u.*sin(poly_loc(:,5)).*sin(poly_loc(:,4)),u.*cos(poly_loc(:,5))];%convert directed motion into cartesean
 

poly_loc(:,1:3)=poly_loc(:,1:3)+U;%iterate translational motion
poly_loc(:,5)=mod(poly_loc(:,5)+omega',pi); %iterate rotation, modulo to constrain between 0 and po
poly_loc(:,4)=mod(poly_loc(:,4)+omega',2*pi); %iterate rotation, modulo to constrain between 0 and po

end