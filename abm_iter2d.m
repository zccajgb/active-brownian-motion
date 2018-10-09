function [poly_loc]= abm_iter2d(PDt,PDr,t,n,poly_loc,u,omega,theta_0); 

randl=[normrnd(0,1,n,3)];%create white noise vectors
titer=[(2*PDt*t)^(0.5)*(randl(:,1:2))];%tranlational brownian iteration of position
riter=[mod((2*PDr*t)^(0.5)*randl(:,3),2*pi)];%rotational browni=an iteration
U=[u.*cos(poly_loc(:,4)+theta_0),u.*sin(poly_loc(:,4)+theta_0)];%convert directed motion into cartesean

poly_loc(:,1:2)=poly_loc(:,1:2)+titer+(U.*t);%iterate translational motion
poly_loc(:,4)=mod(poly_loc(:,4)+riter+(omega'.*t),2*pi); %iterate rotation, modulo to constrain between 0 and po

end