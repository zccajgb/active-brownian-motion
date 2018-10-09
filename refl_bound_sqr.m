function [X,xbndclsn,ybndclsn]= refl_bound_sqr(X,xdist,ydist,xbndclsn,ybndclsn);
%find norm of vector

%index polymersomes outside radius of petri dish   
ind=find(abs(X(:,1))>xdist);
if length(ind)>0
    %reflect polymersome around boundary, modulus is used instead of
    %subtraction incase psome is further out than polymersome radius
      X(ind,1)=sign(X(ind,1)).*(xdist-mod(X(ind,1),xdist));
      

      % count boundary collisions
      xbndclsn=xbndclsn+length(ind); 

end
ind2=find(abs((X(:,2)))>ydist);
if length(ind2)>0
    %reflect polymersome around boundary, modulus is used instead of
    %subtraction incase psome is further out than polymersome radius
      X(ind2,2)=sign(X(ind2,2)).*(ydist-mod(X(ind2,2),ydist));
      

      % count boundary collisions
      ybndclsn=ybndclsn+length(ind2); 
end
end
     
