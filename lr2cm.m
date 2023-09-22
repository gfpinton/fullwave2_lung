%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2020-04-21
% LAST MODIFIED: 2023-09-22
% Launch Fullwave 2 code, zero-pressure approach, easy matlab wrapper
% physical map contains chest wall and lung histology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%addpath /mnt/piave_raid6/oostras/fns
%cd /mnt/piave_raid6/oostras/2023_JASA

cwd=pwd; %current directory to return to that later
nn=1
basedir = ['L12_5_Fd_wmod21_2cm_128' num2str(nn) '/'] %folder where to save txrx data and generated images
eval(['!mkdir -p ' basedir]);

cd(basedir);
input = 'input/' %folder where to save generated acoustical maps
eval(['!mkdir -p ' input]);
cd(cwd);

%for focused sequence:
focal_depth=2e-2 %arbitrarily chosen
L125_spacing_m=1.953e-4 %L12-5_50mm has actual spacing 1.953e-4 [m],
txducer_aperture=64 %number of active txducer elements (each emit)
nelements=192
nevents=nelements-txducer_aperture %walking aperture, sequential txrx events

%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1540;        % speed of sound (m/s)
f0 = 6.25e6;     % frequency [MHz]
omega0 = 2*pi*f0;
lambda = c0/omega0*2*pi
ppw = 12;            % number of points per spatial wavelength
cfl = 0.5;           % Courant-Friedrichs-Levi condition
dX = c0/omega0*2*pi/ppw
dY= c0/omega0*2*pi/ppw
wX = 2.5e-2  % width of simulation field (m), consider angular sensitivity(?)
wY = 3.5e-2         % depth of simulation field (m)
duration = 15.6e-2/c0;  % duration of simulation (s) SHORT TEST!
p0 = 2.5e6;            % confirm!!! pressure in Pa
nX = round(wX/lambda*ppw)  % number of lateral elements
nY = round(wY/lambda*ppw)  % number of depth elements
nT = round(duration*c0/lambda*ppw/cfl)
dT = dX/c0*cfl;
foc=round(focal_depth/dX)      %focal depth in [dX]
fcen=[round(nX/2) foc]         % center of focus in [dX]
beamspacing=round(L125_spacing_m/dX) %spacing in [dX]
fnumber=focal_depth/(txducer_aperture*beamspacing*dX) %asked GP

%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
inmap = zeros(nX,nY); 
inmap(:,1) = ones(nX,1); inmap(:,2) = ones(nX,1); inmap(:,3) = ones(nX,1);
incoords = mapToCoords(inmap);

%%% Generate initial conditions based on input coordinates %%%%%%
ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope
t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icmat=zeros(size(incoords,1),nT);
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
%plot(icvec), hold all
icmat(1:size(incoords,1)/3,:) = focusCoords(fcen(1),fcen(2),incoords(1:size(incoords,1)/3,:),icvec,cfl);
t=t-dT/cfl;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
%plot(icvec)
icmat(size(incoords,1)/3+1:size(incoords,1)/3*2,:)=focusCoords(fcen(1),fcen(2),incoords(size(incoords,1)/3+1:size(incoords,1)/3*2,:),icvec,cfl);
t=t-dT/cfl;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
%plot(icvec), hold off
icmat(size(incoords,1)/3*2+1:size(incoords,1),:)=focusCoords(fcen(1),fcen(2),incoords(size(incoords,1)/3*2+1:size(incoords,1),:),icvec,cfl);
%imagesc(icmat')

%1st
icmat(1:round(nX/2)-round(txducer_aperture*beamspacing/2)-1,:)=0;
icmat(round(nX/2)+round(txducer_aperture*beamspacing/2):nX,:)=0;

%2nd
icmat(nX+1:nX+round(nX/2)-round(txducer_aperture*beamspacing/2)-1,:)=0;
icmat(nX+round(nX/2)+round(txducer_aperture*beamspacing/2):2*nX,:)=0;

%3rd
icmat(2*nX+1:2*nX+round(nX/2)-round(txducer_aperture*beamspacing/2)-1,:)=0;
icmat(2*nX+round(nX/2)+round(txducer_aperture*beamspacing/2):3*nX,:)=0;
%imagesc(icmat')

%%% Generate output coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
outmap = zeros(nX,nY);
outmap(round(nX/2)-round(txducer_aperture*beamspacing/2):round(nX/2)+round(txducer_aperture*beamspacing/2)-1,4) = 1;
outcoords = mapToCoords(outmap);

%%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
materials
rho0 = 1000
A0 = 0.5

load i2365f_etfw1 %body wall

dm=0.33e-3/4; %this is the pixel size in interpd Visual Human slice

%interpolation to elected dX and dY
mat3=interp2easy(cut,dm/dX,0.655*dm/dY,'nearest'); %tissue compression with txducer
mat=zeros(size(mat3,1)+500,size(mat3,2));
mat(1:size(mat3,1),1:size(mat3,2))=mat3;
mat3=mat; clear mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!
crop_depth=round(0.8e-2/dY); %affects the fascia delineation process and depth of lung
%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!
mat3=mat3(crop_depth:crop_depth+nY-1,:)';
xo=dX:dX:dX*size(mat3,1); yo=dX:dX:dX*size(mat3,2);
%imagesc(xo*1e3,yo*1e3,mat3'), axis equal tight, colorbar

cmap=ones(size(mat3))*c0;
rhomap=ones(size(mat3))*rho0;
Amap=ones(size(mat3))*A0;

%assign acoustic properties to each tissue
idc=find(mat3==1 | mat3==4 | mat3==5); %%% CONNECTIVE %%%
cmap(idc)=1613; rhomap(idc)=1120; Amap(idc)=1.57;
idc=find(mat3==2); %%% MUSCLE %%%
cmap(idc)=1547; rhomap(idc)=1050; Amap(idc)=1.09;
idc=find(mat3==3); %%% FAT %%%
cmap(idc)=1478; rhomap(idc)=950; Amap(idc)=0.48;
idc=find(mat3==0); %%% marked FLUID %%%
cmap(idc)=water.c0-40; rhomap(idc)=1000; Amap(idc)=water.alpha;

%smoothing of abd wall maps!!!
gfilt=(5/10)^2*ppw/2;
cmap=imgaussfilt(cmap,gfilt);
rhomap=imgaussfilt(rhomap,gfilt);
Amap=imgaussfilt(Amap,gfilt);

%scatterers generation
nXextend=size(cmap,1);
res_cell = rescell2d(c0,omega0,foc*dX,wX,ncycles,dX,dY);
num_scat = 18;
scat_size = dX;
scat_size_samp = round(scat_size/dX);
scat_lesion = generate_c_scat(1,1,num_scat/res_cell, scat_size_samp, nXextend, nY)-1;
csr=0.035; % scatterer impedance contrast 

air=scat_lesion*0;

%work on the lung histology :
load RL2_m
dm=2.02e-2/814; %size of the original px(?)
lung(1,1)=1; % fix! the single dots of 0 in the bottom corners = 1
lung(1,end)=1;
%figure(1), imagesc(lung), axis equal tight, colorbar
lung=lung(22:end,:); %cut the top border of "nothing"

%interpolation to elected dX and dY
lung=interp2easy(lung,dm/dX,dm/dY,'nearest');
figure(2), imagesc(lung), axis equal tight, colorbar

lungc=ones(size(lung))*c0; % white boards
lungrho=ones(size(lung))*rho0;
lungA=ones(size(lung))*A0;

%ADD VISCERAL PLEURA
vplthickness=0.1e-3; %visceral pleura thickness

for j=1:size(lung,2)
    tissueborder=find(lung(:,j)==0);  % find all the black dots in j-th column
    lung(1:min(tissueborder)-1,j)=-1; % mark the "nothing" above the pleura
    % assign fluid properties to all the shallower dots 
	lungc(1:min(tissueborder)-1,j)=water.c0-40; %marked fluid
    lungrho(1:min(tissueborder)-1,j)=water.rho0;
    lungA(1:min(tissueborder)-1,j)=water.alpha;
	% assign connective tissue properties to the 0.2 mm border
	lungc(min(tissueborder)-round(vplthickness/dY):min(tissueborder)-1,j)=connective.c0;
	lungrho(min(tissueborder)-round(vplthickness/dY):min(tissueborder)-1,j)=connective.rho0;
	lungA(min(tissueborder)-round(vplthickness/dY):min(tissueborder)-1,j)=connective.alpha;
    lung(min(tissueborder)-round(vplthickness/dY):min(tissueborder)-1,j)=-2;
end

originalung=lung(1:round(529*dm/dX),:);
imagesc(originalung), axis equal, axis tight, colorbar

%remove excessive air between pleura and txducer
for j=1:size(lung,2)
    lung(1:find(lung(:,j)==-2,1,'first')-1,j)=-1; % mark the "nothing" above the pleura
    % assign fluid properties to all the shallower dots 
	lungc(find(lung(:,j))==-1,j)=water.c0-40; %marked fluid
    lungrho(find(lung(:,j))==-1,j)=water.rho0;
    lungA(find(lung(:,j))==-1,j)=water.alpha;
end

imagesc(lung), axis equal, axis tight, colorbar

%limit the physical maps by 529 px in depth, no lateral limitations:
lungc=lungc(1:round(529*dm/dX),:);
lungrho=lungrho(1:round(529*dm/dX),:);
lungA=lungA(1:round(529*dm/dX),:);
lung=lung(1:round(529*dm/dX),:);

%%%%%COUNTING AREA ELEMENTS BEFORE MODIFICATIONS%%%%%%%%%%%%%%%%%%
%count area elements of 4 properties: connective tissue, air, soft tissue, water
total_count=size(originalung,1)*size(originalung,2)
pl_fluid_count=size(find(originalung==-1),1)
parench_and_air=total_count-pl_fluid_count
connective_count=size(find(originalung==-2),1)
connective_percent=100*connective_count/parench_and_air
air_count=size(find(originalung==1),1)
air_percent=100*air_count/parench_and_air
soft_tissue_count=size(find(originalung==2),1)
soft_tissue_percent=100*soft_tissue_count/parench_and_air
parench_fluid=size(find(originalung==0),1)
parench_fluid_percent=100*parench_fluid/parench_and_air

%%%%%COUNTING AREA ELEMENTS AFTER MODS%%%%%%%%%%%%%%%%%%
%count area elements of 4 properties: connective tissue, air, soft tissue, water
total_count=size(lung,1)*size(lung,2)
pl_fluid_count=size(find(lung==-1),1)
parench_and_air=total_count-pl_fluid_count
connective_count=size(find(lung==-2),1)
connective_percent=100*connective_count/parench_and_air
air_count=size(find(lung==1),1)
air_percent=100*air_count/parench_and_air
soft_tissue_count=size(find(lung==2),1)
soft_tissue_percent=100*soft_tissue_count/parench_and_air
parench_fluid=size(find(lung==0),1)
parench_fluid_percent=100*parench_fluid/parench_and_air

lungc(find(lung==0))=tissue.c0; %the lung tissue = water/tissue/connective tissue
lungrho(find(lung==0))=tissue.rho0;
lungA(find(lung==0))=tissue.alpha;

lungc(find(lung==1))=340;

%transpose the lung maps:
lungc=lungc';
lungrho=lungrho';
lungA=lungA';

%lung histology positioning:
lat_lung_position=round(size(cmap,1)/2 - size(lungc,1)/2); %lateral
lung_depth=foc-round(1.1e-3/dY); %depth

%finally, add it to the maps!
for i=1:size(lungc,1)
 for j=1:size(lungc,2)
  if lungc(i,j)~=1440
     cmap(lat_lung_position+i-1,lung_depth+j)=lungc(i,j);
	 rhomap(lat_lung_position+i-1,lung_depth+j)=lungrho(i,j);
	 Amap(lat_lung_position+i-1,lung_depth+j)=lungA(i,j);
  end
 end
end

imagesc(cmap'), axis equal tight, colorbar, colormap parula

%lateral subpleural tissue-air interfaces
flayerl=round(0.18e-3/dY);
for i=1:lat_lung_position-1
    pb=find(cmap(i,:)>1500,1,'last');
    cmap(i,pb+flayerl:end)=340;
    cmap(i,pb+flayerl:pb+flayerl+round(vplthickness/dY))=1613;
    rhomap(i,pb+flayerl:pb+flayerl+round(vplthickness/dY))=1120;
    Amap(i,pb+flayerl:pb+flayerl+round(vplthickness/dY))=1.57;
end
imagesc(cmap'), axis equal tight, colorbar, colormap parula
flayerr=round(1.6e-3/dY);
for i=lat_lung_position+size(lungc,1):size(cmap,1)
    pb=find(cmap(i,:)>1500,1,'last');
    cmap(i,pb+flayerr:end)=340;
    cmap(i,pb+flayerr:pb+flayerr+round(vplthickness/dY))=1613;
    rhomap(i,pb+flayerr:pb+flayerr+round(vplthickness/dY))=1120;
    Amap(i,pb+flayerr:pb+flayerr+round(vplthickness/dY))=1.57;
end

%bottom horizontal interface
cmap(:,lung_depth+size(lungc,2)+1:end)=340;

xo=dX:dX:dX*size(cmap,1); yo=dX:dX:dX*size(cmap,2);
imagesc(xo*1e3,yo*1e3,cmap'), axis equal tight, colorbar, colormap parula

air(find(cmap==340))=-1; %air elements

%no scatterers in the air
scat_lesion(find(air==-1))=0;

scat_lesion(find(cmap<water.c0-30 & cmap>350))=0; % no scatterers in the pleural fluid
cmap(find(cmap<water.c0-30 & cmap>350))=water.c0; %pleural fluid = water?
rhomap(find(cmap<water.c0-30 & cmap>350))=water.rho0;
Amap(find(cmap<water.c0-30 & cmap>350))=water.alpha;

%replacement of the air with water in the maps - no propagation in these
%regions, so can be soft tissue properties etc...
cmap(find(air==-1)) = water.c0;
rhomap(find(air==-1)) = water.rho0;
Amap(find(air==-1)) = water.alpha;

rhomap=rhomap-rhomap*rho0*csr; %apply scatterers - rhomap instead of cmap in current version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for event=1:nevents

orig = [round(size(cmap,1)/2)-round(nX/2)-round((event-(nevents+1)/2)*beamspacing)-6 1]

c = cmap(orig:orig+nX-1,1:nY);
rho = rhomap(orig:orig+nX-1,1:nY);
A = Amap(orig:orig+nX-1,1:nY);
air1 = air(orig:orig+nX-1,1:nY);

bovera=c*0-2;
beta=c*0;

xo=-dX*round(size(rho,1)/2):dX:dX*round(size(rho,1)/2); yo=dX:dX:dX*size(rho,2);
imagesc(xo*1e3,yo*1e3,rho'), axis equal tight, colorbar, colormap parula

if event==nevents/2
%%%% show and save acoustical maps
xo=-dX*round(size(c,1)/2):dX:dX*round(size(c,1)/2); yo=dX:dX:dX*size(c,2);
cd([basedir input]);
%c map
imagesc(xo*1e3,yo*1e3,c'), axis equal, axis tight, colorbar, colormap parula,
caxis([1465 1500]),
ax=gca; ax.FontSize=15; ax.FontWeight='Bold';
xlabel('Lateral position [mm]','interpreter','None','fontsize',15,'FontWeight','Bold'),
ylabel('Depth [mm]','interpreter','None','fontsize',15,'FontWeight','Bold')
saveas(gcf, ['cmap_' num2str(event) '.png'], 'png');

%rho map
imagesc(xo*1e3,yo*1e3,rho'), axis equal, axis tight, colorbar, colormap parula
ax=gca; ax.FontSize=15; ax.FontWeight='Bold';
xlabel('Lateral position [mm]','interpreter','None','fontsize',15,'FontWeight','Bold'),
ylabel('Depth [mm]','interpreter','None','fontsize',15,'FontWeight','Bold')
saveas(gcf, ['rhomap_' num2str(event) '.png'], 'png');

%A map
imagesc(xo*1e3,yo*1e3,A'), axis equal, axis tight, colorbar, colormap parula
ax=gca; ax.FontSize=15; ax.FontWeight='Bold';
xlabel('Lateral position [mm]','interpreter','None','fontsize',15,'FontWeight','Bold'),
ylabel('Depth [mm]','interpreter','None','fontsize',15,'FontWeight','Bold')
saveas(gcf, ['Amap_' num2str(event) '.png'], 'png');

%air map
imagesc(xo*1e3,yo*1e3,air1'), axis equal, axis tight, colorbar, colormap gray
ax=gca; ax.FontSize=15; ax.FontWeight='Bold';
xlabel('Lateral position [mm]','interpreter','None','fontsize',15,'FontWeight','Bold'),
ylabel('Depth [mm]','interpreter','None','fontsize',15,'FontWeight','Bold')
saveas(gcf, ['airmap_' num2str(event) '.png'], 'png');

cd(cwd);
end

airinmap = zeros(nX,nY); airinmap(find(air1==-1))=1;
incoordszero = mapToCoords(airinmap);

addpath(cwd);
outdir=[basedir '/txrx_' num2str(event)]; eval(['!mkdir -p ' outdir]); 
eval(['!cp fullwave2_try6_nln_relaxing_pzero ' outdir]);
cd(outdir)
%%%%%%%%%%%%%%%
launch_fullwave2_try6_nln_relaxing4_pzero(c0,omega0,wX,wY,duration,p0,ppw,cfl,c,rho,A,beta,incoords,incoordszero,outcoords,icmat)
cd(cwd);
disp(['txrx ' num2str(event) ' is done'])
end

clear icmat

cd(basedir);
save wmod_L125_2cm_96
cd(cwd);

clear all
disp('launcher is done')