function [poly_loc]= abm_iter3d(PDt,PDr,t,n,poly_loc,u,omega); 

randl=[normrnd(0,1,n,5)];%create white noise vectors
titer=[(2*PDt*t)^(0.5)*(randl(:,1:3))];%tranlational brownian iteration of position
riter=[(2*PDr*t)^(0.5)*randl(:,4:5)];%rotational brownian iteration
U=[u.*sin(poly_loc(:,5)).*cos(poly_loc(:,4)),u.*sin(poly_loc(:,5)).*sin(poly_loc(:,4)),u.*cos(poly_loc(:,5))];%convert directed motion into cartesean

poly_loc(:,1:3)=poly_loc(:,1:3)+titer+(U.*t);%iterate translational motion
poly_loc(:,5)=mod(poly_loc(:,5)+riter(:,1)+(omega'.*t),pi); %iterate rotation, modulo to constrain between 0 and po
poly_loc(:,4)=mod(poly_loc(:,4)+riter(:,2)+(omega'.*t),2*pi); %iterate rotation, modulo to constrain between 0 and po

end