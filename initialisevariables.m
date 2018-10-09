function [k,T,eta,Dtg,petriR,polyR,PDt,PDr,bndclsn,pcolsn,rad,azi,eli,msd]=initialisevariables(M,loops)

k=1.38064852e-23;%boltzmann constant
% k=1;
T=298.15%temperature
eta=8.90e-4%viscosity
 Dtg=7e-10;%Longsworth, L. G. 1955. Diffusion in liquids and the Stokes-Einstein relation, p. 225-247. In T. Shedlovsky (ed.), Electrochemistry in biology and medicine. John Wiley & Sons, Inc., New York, N.Y.
% Dtg=0.7;

petriR=0.05; %radius of petri dish
polyR=50e-9; %radius of polymersome 
zax=1e-3;%z axis of polymersome
% poly_loc=psome_init(M,petriR,zax); %randomly initialises polymersomes
 %translational diffusion coeff for the psome
% PDt=0.22e-6;
% PDr=0.16; 
 PDt=(k.*T)./(6.*pi.*eta.*polyR); %translational diffusion coeff of poly
PDr=(k.*T)./(8.*pi.*eta.*polyR.^3);%rotational diffusion coeff for the psome
bndclsn=0;%counter for number of colisions with boundary
pcolsn=0;%counter for number of polymersome colisions
%radial,azimuthal and elivation tracker for histogram
rad=zeros(M,loops);
azi=zeros(M,loops);
eli=zeros(M,loops);
msd=0;
end