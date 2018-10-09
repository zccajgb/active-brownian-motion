function [k,T,eta,Dtg,petriR,polyR,PDt,PDr,bndclsn,pcolsn,rad,azi,eli,msd]=initialisevariablesGs(M,loops)

k=1.38064852*10-5;%boltzmann constant
% k=1;
T=298.15%temperature
eta=8.90e5%viscosity
% Dtg=.7e-10;%Longsworth, L. G. 1955. Diffusion in liquids and the Stokes-Einstein relation, p. 225-247. In T. Shedlovsky (ed.), Electrochemistry in biology and medicine. John Wiley & Sons, Inc., New York, N.Y.
Dtg=0.7;
%Dt=(k*T)/(6*pi*eta*r); %translational diffusion coeff of glucose
petriR=0.05; %radius of petri dish
polyR=0.05e-6; %radius of polymersome 
zax=1e-3;%z axis of polymersome
% poly_loc=psome_init(M,petriR,zax); %randomly initialises polymersomes
PDt=220 %translational diffusion coeff for the psome
PDr=1.6e8%rotational diffusion coeff for the psome
bndclsn=0;%counter for number of colisions with boundary
pcolsn=0;%counter for number of polymersome colisions
%radial,azimuthal and elivation tracker for histogram
rad=zeros(M,loops);
azi=zeros(M,loops);
eli=zeros(M,loops);
msd=0;
end
