function [X,bndclsn]= refl_bound3d(X,R,bndclsn);
%find norm of vector
nrml=sqrt(sum(abs(X(:,1:3)).^2,2));
%index polymersomes outside radius of petri dish   
ind=find((nrml)>R);
if length(ind)>0
    %reflect polymersome around boundary, modulus is used instead of
    %subtraction incase psome is further out than polymersome radius
      r=R-mod(nrml(ind),R);
      %find polar and azimuthal angles
      [az,el,~]=cart2sph(X(ind,1),X(ind,2),X(ind,3));
      %convert to spherical coordintates
      [xi,eta,zeta]=sph2cart(az,el,r);
      %update polymersome location
      X(ind,1:3)=[xi,eta,zeta];
      % count boundary collisions
      bndclsn=bndclsn+length(ind); 
end
     