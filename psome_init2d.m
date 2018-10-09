function [poly_locs]=psome_init2d(M,R) 
% phi=2*pi*rand(M,1); %(uniform) randomly generates an azmiuth between 0-2pi
% theta=pi*rand(M,1); %randomly generate polar between 0-pi
% rho=rand(M,1)*R; %(uniform) randomly generates a radius between 0-R
x=2*R*rand(M,1)-R;
y=2*R*rand(M,1)-R;
z=zeros(M,1);
% z=2*R*rand(M,1)-R;
% z=2*zax*rand(M,1)-zax;
 inclin=zeros(M,1);

azimuth=2*pi*rand(M,1);

% [x,y,z]=pol2cart(phi,theta,rho); %converts to cartesean

poly_locs=[x,y,z,azimuth,inclin]; %converts to complex cartesean

end