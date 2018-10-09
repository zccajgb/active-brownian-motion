function [X,bndclsn]= refl_bound_cyl(X,R,bndclsn,zax);
%find norm of vector
nrml=sqrt(sum(abs(X(:,1:2)).^2,2));

%index polymersomes outside radius of petri dish   
ind=find((nrml)>R);
ind2=find(abs(X(:,3))>zax);
if length(ind)>0
    %reflect polymersome around boundary, modulus is used instead of
    %subtraction incase psome is further out than polymersome radius
      r=R-mod(nrml(ind),R);
      %find polar and azimuthal angles
      [theta,~]=cart2pol(X(ind,1),X(ind,2));
      %convert to spherical coordintates
      [xi,eta]=pol2cart(theta,r);
      %update polymersome location
      X(ind,1:3)=[xi,eta,X(ind,3)];
      % count boundary collisions
      bndclsn=bndclsn+length(ind); 
end

if length(ind2)>0
    %reflect polymersome around boundary, modulus is used instead of
    %subtraction incase psome is further out than polymersome radius
      ind3=find(X(ind2,3)>0);
      Z3=zax-mod(X(ind3,3),zax);
      ind4=find(X(ind2,3)<0);
      Z4=-zax-mod(X(ind4,3),zax);
      %find polar and azimuthal angles
      [theta3,r3]=cart2pol(X(ind3,1),X(ind3,2));
       [theta4,r4]=cart2pol(X(ind4,1),X(ind4,2));
      %convert to spherical coordintates
      [xi3,eta3]=pol2cart(theta3,r3);
      [xi4,eta4]=pol2cart(theta4,r4);
      %update polymersome location
      X(ind3,1:3)=[xi3,eta3,Z3];
      X(ind4,1:3)=[xi4,eta4,Z4];
      % count boundary collisions
      bndclsn=bndclsn+length(ind);
end
end