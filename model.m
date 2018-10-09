beep off
close all
%choose mode and number of polyersomes
%mode determines method of calculating velocities
clear
disp('Select mode:')
disp('1 = constant mu, non-constant G, ASA')
disp('2 = non-constant mu,non-constant G, ASA for u0, ASA for u1')
disp('0 = test')

mode=input('');
M=input('Number of polymersomes: ');
makeanimation=input('Make animation: 0 for no, 1 for yes: ');
geom=input('Petri Dish (0) or Rectangle (1): ');
%time step in 1*10e^9 s
t=0.005
ploc=[];
% %repeatable random numbers for debugging
% rng('default');
% rng(1);
%random random numbers for simulation
rng('shuffle');
%number of loo
loops=(200);
tau0=t*loops;
% a=-10:1:10;
% mup=[10.^a];
mup=-1e-8
tic
for index17=1:length(mup)
    tic
    mupeo=mup(index17)
    if mode==1
    mupmpc=mupeo;
    elseif mode==2
    mupmpc=0;
    end
%Initialise variables (in a seperate function for tidyness)
[k,T,eta,Dtg,petriR,polyR,PDt,PDr,bndclsn,pcolsn,rad,azi,eli,msd]=initialisevariables(M,loops);
[kpmpc,kpeo,theta_0,phi_0]=initialisevariables2;
xbndclsn=0;
ybndclsn=0;
poly_loc=psome_init2d(M,petriR); %randomly initialises polymersomes]
if geom==1
xdist=0.0285;
ydist=0.0092;
initposit=[-0.021*ones(M,1),zeros(M,1)];
% initposit=([zeros(M,2)])
poly_loc(:,1:2)=initposit;
end

%create file for animation
H(loops) = struct('cdata',[],'colormap',[]);


for i=1:(loops)
i
timeiter=i;
%different modes
for klklk=1
if mupeo==0
    u=0;
    omega=0;
else
if mode==0 %test, arbitrary values
    u=0;
    omega=0;
elseif mode==1 %const mu, nc gluc, asa1
    omega=0;
    Gpoly=0;
    
    kpmpc=0;
    u=cmcgasa(poly_loc,polyR,Dtg,kpeo,kpmpc,mupeo,mupmpc,timeiter,t,phi_0,theta_0,M,Gpoly,tau0);
elseif mode==2 %nc mu, mc gluc, asa
    deltamu=mupeo-mupmpc;
    Gpoly=0;
    omega=rotngasa(poly_loc,polyR,Dtg,kpeo,kpmpc,deltamu,timeiter,t,phi_0,theta_0,M,Gpoly,petriR);
    u0=cmngasa(poly_loc,polyR,Dtg,kpeo,kpmpc,mupeo,mupmpc,timeiter,t,phi_0,theta_0,M,Gpoly,tau0);
    u1=nmngasa(poly_loc,polyR,Dtg,kpeo,kpmpc,deltamu,timeiter,t,phi_0,theta_0,M,Gpoly);
    u=u0+u1';
end
end
end
%iterates movement of polymersomes
x=poly_loc(:,1);
y=poly_loc(:,2);
ploc{i}=poly_loc;
% z=poly_loc(:,3);
[poly_loc]= abm_iter2d(PDt,PDr,t,length(poly_loc(:,1)),poly_loc,u,omega,theta_0);
%imposes reflective boundary conditions on polyermesomes
if geom==0
[poly_loc,bndclsn]=refl_bound2d(poly_loc,petriR,bndclsn);
elseif geom==1
[poly_loc,xbndclsn,ybndclsn]= refl_bound_sqr(poly_loc,xdist,ydist,xbndclsn,ybndclsn);
end
%account for polymersome collisions
% [poly_loc]=pp_colsn(poly_loc,polyR, pcolsn); 
%calculate msd of step
samptime(i,index17)=i.*t;
% msdloop=mean((poly_loc(:,1)-x).^2+(poly_loc(:,2)-y).^2+(poly_loc(:,3)-z).^2);
% displacement=[(poly_loc(:,1)-x),(poly_loc(:,2)-y),(poly_loc(:,3)-z)];
msdloop=mean((poly_loc(:,1)-x).^2+(poly_loc(:,2)-y).^2)./t;
medmsdloop=median((poly_loc(:,1)-x).^2+(poly_loc(:,2)-y).^2)./t;
displacement=[(poly_loc(:,1)-x),(poly_loc(:,2)-y)];
vel(:,:)=displacement./samptime(i,index17);
% driftvel1(:)=sqrt(vel(:,1).^2+vel(:,2).^2+vel(:,3).^2);
driftvel1(:)=sqrt(vel(:,1).^2+vel(:,2).^2);
driftvel(i,index17)=mean(driftvel1);
meddriftvel(i,index17)=mean(driftvel1);
% if i==1
%    msdloop=0
% end
msd=msd+msdloop;
%create data for histograms
[azi(:,i),rad(:,i)]=cart2pol(poly_loc(:,1),poly_loc(:,2));
for plottingind=1
%plot histograms on every iteration
% subplot(2,2,1)
% histogram(rad(:,i))
% xlabel('R')
% subplot(2,2,2)
% histogram(azi(:,i))
% xlabel('azimuth')
% subplot(2,2,3)
% histogram(eli(:,i))
% xlabel('elevation')
% subplot(2,2,4)
% cla
% viscircles([0,0],petriR,'color','black','linewidth',0.5) %plot boundary of petri dish
%    hold on
% scatter((poly_loc(:,1)),(poly_loc(:,2))) %plot psomes
% axis([-petriR,petriR,-petriR,petriR]);
% xlabel('x')
% ylabel('y')
%% animation
if makeanimation==1
    
%2d plot of x and y
% subplot(2,2,1)
if geom==0
cla
%plot boundary of petri dish
viscircles([0,0],petriR,'color','black','linewidth',0.5)
   hold on
% plot polymersomes
scatter(initposit(1,1),initposit(1,2),36,'red')
scatter((poly_loc(:,1)),(poly_loc(:,2)),36,'blue') 
axis([-petriR,petriR,-petriR,petriR]);
xlabel('x')
ylabel('y')
elseif geom==1
cla
%plot boundary of petri dish
%rectangle('Position',[-xdist,-ydist,2*xdist,2*ydist],'EdgeColor','black','linewidth',0.5)
   hold on
% plot polymersomes
%scatter(initposit(1,1),initposit(1,2),36,'red','+')
scatter(0,0,36,'red','+')
%scatter(0,0,36,'black','x')
scatter((poly_loc(:,1)-initposit(1,1)),(poly_loc(:,2)-initposit(1,2)),36,'blue')
% legend('Initial Polyersome Position','Initial Glucose Position','Polymersomes')
legend('Initial Polyersome Position','Polymersomes')
axis([-2e-5,2e-5,-2e-5,2e-5]);
xlabel('x')
ylabel('y')
end

%2d plot of x and z
% subplot(2,2,2)
%  cla
% viscircles([0,0],petriR,'color','black','linewidth',0.5) %plot boundary of petri dish
%    hold on
% scatter((poly_loc(:,1)),(poly_loc(:,3))) %plot psomes
% axis([-petriR,petriR,-petriR,petriR]);
str=sprintf('%d',i);
title(str)
% xlabel('x')
% ylabel('z')

%2d plot of y and z
% subplot(2,2,3)
% cla
% viscircles([0,0],petriR,'color','black','linewidth',0.5) %plot boundary of petri dish
%    hold on
% scatter((poly_loc(:,2)),(poly_loc(:,3))) %plot psomes
% axis([-petriR,petriR,-petriR,petriR]);
% xlabel('y')
% ylabel('z')
% 
% %3d plot
% subplot(2,2,4)
% cla
% scatter3(poly_loc(:,1),poly_loc(:,2),poly_loc(:,3))
% axis([-petriR,petriR,-petriR,petriR]);
% view(3)
% forces plots to be made on each loop rather than waiting to the end
drawnow
% take snapshot of plot in order to make animation
H(i)=getframe;

cla

% for larger loops this was used to take data every minute to have smaller
% data set
end
end  

% if mod(i,10)==0

 msdout(i,index17)=msdloop;  

  tu(:,i,index17)=u(:);
% end
end
%create gif file from movie
if makeanimation==1
for i=1:loops
gifmaker(H(i),i,'animation.gif')
end
end
toc
actualmsd(index17)=sum(msdout(:,index17))./(t*i)
% if abs(actualmsd(index17)-100e-6)<10e-6
%     mup(index17)
%     break
% end
toc
end
close all
% driftvelocity(:,index17)=msd./(samptime(:,index17).*1e9);
%plot histograms
% subplot(1,2,1)
histogram(rad(:))
xlabel('R')
% subplot(1,2,2)
% histogram(azi(:))
% xlabel('azimuth')
% subplot(2,2,3)
% histogram(eli(:))
% xlabel('elevation')
% subplot(2,2,4)
if makeanimation==1
% savefig('histogram')
save('H.mat','H')
end
% r=20e-6;
% d = r*2;
% px = 0.04-r;
% py = 0-r;
% h = rectangle('Position',[px py d d],'Curvature',[1,1]);
% pos = [0.0002 0.0004 0.002 0.002];
% rectangle('Position',pos,'Curvature',[1 1])

% movie(H,1,1000)
% for m=1:k
% viscircles([0,0],R,'color','black','linewidth',0.5) %plot boundary of petri dish
%    hold on
%    scatter(real(pso{m}(:,1)),imag(pso{m}(:,1))) %plot psomes
% 
% ind=find(gluc{m}(:,3)==0)
%    scatter(real(gluc{m}(ind,1)),imag(gluc{m}(ind,1)),'.') %plot glucose
% ind=find(X(:,3)==1)
%    scatter(real(gluc{m}(ind,1)),imag(gluc{m}(ind,1)),'.') %plot d-lactone
%    pause
%  hold off
% end

    %axis([-0.1,0.1,-0.1,0.1])