%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2020-04-21
% LAST MODIFIED: 2021-04-21
% Launch Fullwave 2 code, zero-pressure approach, easy matlab wrapper
% hybrid physical map contains abdominal wall and synthetic 2x5 mm cluster of staggered circular air bubbles of radius=dX (which is minimal possible for used functions),
% no tissue walls between them:
% bubble centers are tied, lateral distance between centers (306 microm);
% axial distance between horizontal levels of the centers (166 microm);
% the lesion is surrounded by flat tissue-air interface (laterally and at the bottom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
addpath /pine/scr/o/l/oleksii
%addpath /mnt/intra_raid6/oostras
cd /pine/scr/o/l/oleksii
%cd /mnt/intra_raid6/oostras
basedir = 'syn_ntdX_128/' %folder where to save all txrx data
eval(['!mkdir -p ' basedir]);

%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1540;         % speed of sound (m/s)
omega0 = 2*pi*5e6; % center radian frequency of transmitted wave
wX = 3e-2;         % width of simulation field (m)
wY = 2.7e-2;         % depth of simulation field (m)
duration = wY*5.8/c0;  % duration of simulation (s)
p0 = 1e5; % pressure in Pa

%%% Advanced variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ppw = 12;           % number of points per spatial wavelength
cfl = 0.5;         % Courant-Friedrichs-Lewy condition

%%% Grid size calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = c0/omega0*2*pi;
nX = round(wX/lambda*ppw);  % number of lateral elements
nY = round(wY/lambda*ppw);  % number of depth elements
nT = round(duration*c0/lambda*ppw/cfl);
dX = c0/omega0*2*pi/ppw
dY= c0/omega0*2*pi/ppw
dT = dX/c0*cfl;
foc=round(1.9e-2/dX);    %focal depth
fcen=[round(nX/2) foc];  % center of focus

%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
inmap = zeros(nX,nY); 
inmap(:,1) = ones(nX,1); inmap(:,2) = ones(nX,1); inmap(:,3) = ones(nX,1);
incoords = mapToCoords(inmap);

%%% Generate initial conditions based on input coordinates %%%%%%
ncycles = 2; % number of cycles in pulse
dur = 2;     % exponential drop-off of envelope
t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icmat=zeros(size(incoords,1),nT);
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
icmat(1:size(incoords,1)/3,:) = focusCoords(fcen(1),fcen(2),incoords(1:size(incoords,1)/3,:),icvec,cfl);
t=t-dT/cfl;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
icmat(size(incoords,1)/3+1:size(incoords,1)/3*2,:)=focusCoords(fcen(1),fcen(2),incoords(size(incoords,1)/3+1:size(incoords,1)/3*2,:),icvec,cfl);
t=t-dT/cfl;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
icmat(size(incoords,1)/3*2+1:size(incoords,1),:)=focusCoords(fcen(1),fcen(2),incoords(size(incoords,1)/3*2+1:size(incoords,1),:),icvec,cfl);

%%% Generate output coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
outmap = zeros(nX,nY);  outmap(:,4) = ones(nX,1);
outcoords = mapToCoords(outmap);

%%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
materials
rho0=1000;
A0 = 3/(20/log(10));
N0 = -tissue.beta/(rho0*c0^4);
fnumber=foc/nX;
beamwidth=round(lambda*fnumber/dX);
beamspacing=beamwidth/2;
nlines=128;
[m.c m.rho m.A m.N] = img2fieldFlatten('r102gh.tif',dX,dY,c0,rho0);

%abdominal wall added to the physical maps
cext=ones(size(m.c,1), nY)*c0;
cext(:, 1:size(m.c, 2))=m.c;
rhoext=ones(size(m.rho,1), nY)*rho0;
rhoext(:, 1:size(m.rho, 2))=m.rho;
Aext=ones(size(m.A,1), nY)*A0;
Aext(:, 1:size(m.A, 2))=m.A;
Next=ones(size(m.N,1), nY)*N0;
Next(:, 1:size(m.N, 2))=m.N;

%add fascia and parietal pleura
pplthickness=4e-4;      %(0.4 mm thick layer)
cext(:, round((1.9e-2-pplthickness)/dY):round((1.9e-2)/dY)) = connective.c0;
rhoext(:, round((1.9e-2-pplthickness)/dY):round((1.9e-2)/dY)) = connective.rho0;
Aext(:, round((1.9e-2-pplthickness)/dY):round((1.9e-2)/dY)) = connective.alpha;

%mark the cmap deeper than parietal pleura
cext(:, round((1.9e-2)/dY)+1:end) = water.c0+1;

%scatterers generation
nXextend=size(m.c,1);
res_cell = rescell2d(c0,omega0,foc*dX,wX,ncycles,dX,dY);
num_scat = 20;
scat_size = dY;
scat_size_samp = round(scat_size/dX);
cscat_lesion = generate_c_scat(1,1,num_scat/res_cell, scat_size_samp, nXextend, nY)-1;
csr=0.05; % scatterer impedance contrast 

air=cscat_lesion*0; %matrix for the air 

%where to place the structure?
for n=1:nlines
range(n)=round(3.5e-2/dX-(n-(nlines+1)/2)*beamspacing);
end
ledge=min(range); redge=(max(range));               % left and right edges of the B-mode image in the future
sposition = round(ledge+nX/2+(redge-ledge)/2);      % lateral position of the base center (X) - center of the bmode image now

%add a bell of bare connective tissue surrounded by the air
basewidth=2e-3; % the widest bell size (X)
sdepth=5e-3;    % length/depth of the bell (Y)

ncentral = round(sdepth/dY);   % number of elements (Y) in the central vertical column
nbase = round(basewidth/dX/2); % number of elements (X) in a base row on each side from the central column

for j=1:ncentral               % depth counter (Y)
  for i=-nbase:nbase           %lateral counter (X)
    if    nbase-abs(i)>exp(4*j/ncentral)
    cext(sposition+i, round(1.9e-2/dY)+j) = connective.c0;
	end
  end
  %find tissue border in j-th row    
    tissueborder=find(cext(:,round(1.9e-2/dY)+j)==connective.c0);
	lborder=min(tissueborder); rborder=max(tissueborder);
    %add 8px of air laterally to the tissue border
    air(lborder-8:lborder, round(1.9e-2/dY)+j)=-1;
	air(rborder:rborder+8, round(1.9e-2/dY)+j)=-1;
end

cscat_lesion(find(cext==water.c0+1))=0;

%horizontal straight tissue-air interfaces
air(1:sposition-nbase, round(1.9e-2/dY)+1:round(1.9e-2/dY)+8)=-1;   %left
air(sposition+nbase:end, round(1.9e-2/dY)+1:round(1.9e-2/dY)+8)=-1; %right

%augmentation of the curve extremum (air layer thickness at least 8px along the whole border)
air(sposition-8:sposition+8, max(find(cext(sposition,:)==connective.c0))+1:max(find(cext(sposition,:)==connective.c0))+8)=-1;

%do NOT add large water oval bubbles
%add small oval air bubbles in the center of the water bubbles
wbubbleradius=14e-5;             %lateral (X) radius of the bubbles, as exudate-filled part of infiltrated alveoli
waxialradius=wbubbleradius/2;   %axial (Y) radius of the bubbles =270um diam
wbubbledistance=dX;           %distance between water bubbles (borders, not centers), dX multiplicity!
bubbleradius=dX;              %radius of the air bubbles, as aerated parts of infiltrated alveoli

for j=round((1.9e-2+waxialradius)/dY):round(2*waxialradius/dY)+1:round((1.9e-2-waxialradius)/dY)+ncentral-1  
  for i=sposition-nbase:round((wbubbledistance+2*wbubbleradius)/dX):sposition+nbase
     thumbup=0;
	 if (mod(abs(j-round((1.9e-2+waxialradius)/dY))/(round(2*waxialradius/dY)+1), 2)~=1) & (cext(i, j)==connective.c0)
        ida=ovalIdx(size(cscat_lesion),[i j],round(bubbleradius/dX),round(bubbleradius/dX/2));
        thumbup=thumbup+1;
     end		
     if (mod(abs(j-round((1.9e-2+waxialradius)/dY))/(round(2*waxialradius/dY)+1), 2)==1) & (cext(i+round((wbubbledistance+2*wbubbleradius)/dX/2), j)==connective.c0)
        ida=ovalIdx(size(cscat_lesion),[i+round((wbubbledistance+2*wbubbleradius)/dX/2) j],round(bubbleradius/dX),round(bubbleradius/dX/2));
		thumbup=thumbup+1;
     end
	    if thumbup>0
        air(ida)=-1;
     end
   end
end

%no scatterers in the air
cscat_lesion(find(air==-1))=0;

%replace the map deeper than parietal pleura with water
cext(:, round((1.9e-2)/dY)+1:end) = water.c0;
rhoext(:, round((1.9e-2)/dY)+1:end) = water.rho0;
Aext(:, round((1.9e-2)/dY)+1:end) = water.alpha;

%smoothing of abd wall maps
gfilt=(5/10)^2*ppw/2;
cext=imgaussfilt(cext,gfilt);
rhoext=imgaussfilt(rhoext,gfilt);
Aext=imgaussfilt(Aext,gfilt);

for n=1:nlines
orig = [round(3.5e-2/dX-(n-(nlines+1)/2)*beamspacing) 1]
  c = chopField(cext,c0,orig,nX,nY);
  rho = chopField(rhoext,rho0,orig,nX,nY);
  A = chopField(Aext,A0,orig,nX,nY);
  N = chopField(Next,N0,orig,nX,nY);
   
  cs = cscat_lesion(orig(1):orig(1)+nX-1,:);
  air2 = air(orig(1):orig(1)+nX-1,:);

%remove the air at the edges!
  air2(1:round(size(air2,1)/20),:)=0;
  air2(size(air2,1)-round(size(air2,1)/20):size(air2,1),:)=0;
  air2(:,1:round(size(air2,2)/20))=0;
  air2(:,size(air2,2)-round(size(air2,2)/20):size(air2,2))=0;

  bovera=N*0-2;
  beta=c*0;

  c=c-cs*c0*csr;
  
  airinmap = zeros(nX,nY); airinmap(find(air2==-1))=1;
  incoordszero = mapToCoords(airinmap);
  
if n==round(nlines/2)
prevdir=pwd;
cd(basedir);
  % show and save c map
latm=-nX*dX/2:dX:nX*dX/2; depm=dY:dY:nY*dY;
imagesc(latm*1e3,depm*1e3,c'), axis equal, axis tight, colorbar, colormap parula
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title('cmap, 5 MHz, 12 ppw', 'Interpreter', 'None', 'FontSize', 10, 'FontWeight', 'normal')
saveas(gcf, 'cmap_n.png', 'png');

%cmap cropped to b-mode borders
cbmap=c(round(nX/2-n*beamspacing):round(nX/2+n*beamspacing), 3:end);
latbm=-size(cbmap, 1)*dX/2:dX:size(cbmap, 1)*dX/2; depbm=dY*3:dY:nY*dY;
imagesc(latbm*1e3,depbm*1e3,cbmap'), axis equal, axis tight, colorbar, colormap parula
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title('cmap, visible part, 5 MHz, 12 ppw', 'Interpreter', 'None', 'FontSize', 10, 'FontWeight', 'normal')
saveas(gcf, 'cmap_visible.png', 'png');

%zoomed cmap
czmap=c(round(nX/2)-round(2e-3/dX):round(nX/2)+round(2e-3/dX), round(1.85e-2/dY):round(2.4e-2/dY));
latzm=-size(czmap, 1)*dX/2:dX:size(czmap, 1)*dX/2; depzm=1.85e-2:dY:2.4e-2;
imagesc(latzm*1e3,depzm*1e3,czmap'), axis equal, axis tight, colorbar, colormap parula
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title('zoomed cmap, 5 MHz, 12 ppw', 'Interpreter', 'None', 'FontSize', 10, 'FontWeight', 'normal')
saveas(gcf, 'cmap_zoomed.png', 'png');

 % show and save rho map
imagesc(latm*1e3,depm*1e3,rho'), axis equal, axis tight, colorbar, colormap parula
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title('rhomap, 5 MHz, 12 ppw', 'Interpreter', 'None', 'FontSize', 10, 'FontWeight', 'normal')
saveas(gcf, 'rhomap_n.png', 'png');

%rhomap cropped to b-mode borders
rhobmap=rho(round(nX/2-n*beamspacing):round(nX/2+n*beamspacing), 3:end);
imagesc(latbm*1e3,depbm*1e3,rhobmap'), axis equal, axis tight, colorbar, colormap parula
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title('rhomap, visible part, 5 MHz, 12 ppw', 'Interpreter', 'None', 'FontSize', 10, 'FontWeight', 'normal')
saveas(gcf, 'rhomap_visible.png', 'png');

 % show and save air map
imagesc(latm*1e3,depm*1e3,air2'), axis equal, axis tight, colorbar, colormap parula
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title('air map, 5 MHz, 12 ppw', 'Interpreter', 'None', 'FontSize', 10, 'FontWeight', 'normal')
saveas(gcf, 'airmap_n.png', 'png');

%air map cropped to b-mode borders
airbmap=air2(round(nX/2-n*beamspacing):round(nX/2+n*beamspacing), 3:end);
imagesc(latbm*1e3,depbm*1e3,airbmap'), axis equal, axis tight, colorbar, colormap parula
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title('air map, visible part, 5 MHz, 12 ppw', 'Interpreter', 'None', 'FontSize', 10, 'FontWeight', 'normal')
saveas(gcf, 'airmap_visible.png', 'png');

%zoomed airmap
airzmap=air2(round(nX/2)-round(2e-3/dX):round(nX/2)+round(2e-3/dX), round(1.85e-2/dY):round(2.4e-2/dY));
imagesc(latzm*1e3,depzm*1e3,airzmap'), axis equal, axis tight, colorbar, colormap parula
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title('zoomed air map, 5 MHz, 12 ppw', 'Interpreter', 'None', 'FontSize', 10, 'FontWeight', 'normal')
saveas(gcf, 'airmap_zoomed.png', 'png');

 % show and save scatteres map
imagesc(latm*1e3,depm*1e3,cs'), axis equal, axis tight, colorbar, colormap parula
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title('scatterers map, 5 MHz, 12 ppw', 'Interpreter', 'None', 'FontSize', 10, 'FontWeight', 'normal')
saveas(gcf, 'scatmap_n.png', 'png');

%scatterers map cropped to b-mode borders
scatbmap=cs(round(nX/2-n*beamspacing):round(nX/2+n*beamspacing), 3:end);
imagesc(latbm*1e3,depbm*1e3,scatbmap'), axis equal, axis tight, colorbar, colormap parula
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title('scatterers map, visible part, 5 MHz, 12 ppw', 'Interpreter', 'None', 'FontSize', 10, 'FontWeight', 'normal')
saveas(gcf, 'scatmap_visible.png', 'png');
cd(prevdir);
end

  cwd=pwd; addpath(cwd);
  outdir=[basedir '/txrx_' num2str(n)]; eval(['!mkdir -p ' outdir]); 
  eval(['!cp fullwave2_try6_nln_relaxing_pzero ' outdir]);
  cd(outdir)
  launch_fullwave2_try6_nln_relaxing4_pzero(c0,omega0,wX,wY,duration,p0,ppw,cfl,c,rho,A,beta,incoords,incoordszero,outcoords,icmat)
  cd(cwd);
n
disp('is done')
end

clear icmat

cd(basedir);

save syn_ntdX_128

cd(prevdir);

clear all
disp('launcher is done')