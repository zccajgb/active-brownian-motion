function [poly_loc]= pp_colsn(poly_loc, polyR, pcolsn)
for N=1:length(poly_loc)  
p2=poly_loc-poly_loc(N);
p2(N)=10; %arbritarily large number
ind=find(abs(p2)<2*polyR);
if length(ind)>0
       Rvec=ones(length(ind),1)*polyR;
      r=Rvec+0.5*abs(p2(ind));
      theta=angle(p2(ind));
      [rpart,ipart]=pol2cart(theta,polyR);
      poly_locs(ind)=rpart+ipart*i;
      pcolsn=pcolsn+length(ind);
end
end
end
