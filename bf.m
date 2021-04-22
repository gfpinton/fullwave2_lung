%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2019-11-21
% LAST MODIFIED: 2020-08-20
% Basic transducer data read and beamform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
addpath /pine/scr/o/l/oleksii
%addpath /mnt/intra_raid6/oostras
cd /pine/scr/o/l/oleksii
%cd /mnt/intra_raid6/oostras

basedir = 'syn_ntdX_128/'
cd(basedir)
load syn_ntdX_128

cd /pine/scr/o/l/oleksii
%cd /mnt/intra_raid6/oostras

%% This is a raw single TX-RX event tt=1:128
tt=1; % 
outdir=[basedir 'txrx_' num2str(tt)]; 
nRun=sizeOfFile([outdir '/genout.dat'])/4/size(outcoords,1)
pxducer=readGenoutSlice([outdir '/genout.dat'],0:nRun-2,size(outcoords,1));
%figure(1), imagesc(powcompress(pxducer,1/4))

%% The transducer coordinates 
xducercoords=outcoords;

%addpath /celerina/gfp/mfs/fullwave_tissue_calibration/
pxducer2=gen_txrx_pxducer(basedir,size(outcoords,1),1:size(outcoords,1),nT,nlines,1:size(outcoords,1),ones(nT,1));
%figure(2), imagesc(powcompress(pxducer2(:,:,1),1/4))

%%% GENERATE IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deps = dY*3:lambda/8:0.9*nT*dT*c0/2;
lats = 0;
xducercoords = outcoords;
bm=zeros(length(lats),length(deps),nlines);

%% Generate beamforming table %%
idt0=0; % WITHOUT IDT0! %
idps_mat=zeros(length(lats),length(deps),nT+1);
for ii=1:length(lats)
  lat=lats(ii);
  for jj=1:length(deps)
    dep=deps(jj);
    fcen=round([lat/dY+mean(xducercoords(:,1)) dep/dY ]);
    idx=find(abs(xducercoords(:,1)-fcen(1))<=fcen(2)/fnumber);
    dd=focusProfile(fcen,xducercoords(idx,:),dT/dY*c0);
    idt=idt0+round(2*dep/double(c0)/(dT));
    idp=double((size(pxducer2,1)*(idx-1))+double(idt)+dd);
    idps_mat(ii,jj,1)=length(idp);
    idps_mat(ii,jj,2:length(idp)+1)=(idp);
  end
end

%% delay and sum %%
ii=1;
pxducer3=pxducer2(:,:,1);
%px=pxducer3(:,round(end/2)); %'end/2' did not work
px=pxducer3(:,round(size(pxducer3,2)/2));
[val idt0]=max(abs(hilbert(px)))
for n=1:nlines
  pxducer3=pxducer2(:,:,nlines-n+1);
  for jj=1:length(deps)
    bm(ii,jj,n)=sum(pxducer3(round(idps_mat(ii,jj,2:idps_mat(ii,jj,1)))+idt0));
  end
end

%to enable saving the workspace
clear idps_mat

cd(basedir);
save syn_ntdX_128_bmode

%%RENDERING%%
%show and save B-mode image
n=1:nlines; bws=((n-(nlines+1)/2)*beamspacing)*dY;
imagesc(bws*1e3,deps*1e3,dbah(squeeze(bm)),[-60 0])
ax=gca;
ax.FontSize=15;
ax.FontWeight='Bold';
colormap gray, colorbar
xlabel('Lateral position [mm]','interpreter','None','fontsize',15,'FontWeight','Bold'),
ylabel('Depth [mm]','interpreter','None','fontsize',15,'FontWeight','Bold')
axis equal, axis tight
%title('bell shaped lesion filled with air bubbles in plain water, 5 MHz, 12 ppw', 'Interpreter', 'None', 'FontSize', 10, 'FontWeight', 'normal')
saveas(gcf, 'syn_ntdX_128_bmode.svg', 'svg');
saveas(gcf, 'syn_ntdX_128_bmode.png', 'png');

% BASIC TGC %
for n=1:nlines
 for j=1:size(bm,2)
   bm_tgc(:,j,n)=bm(:,j,n)*exp(2*100*deps(j)/8.686); %doubled TGC  
 end
end

%%RENDERING WITH APPLIED TGC%%
%show and save B-mode image
n=1:nlines; bws=((n-(nlines+1)/2)*beamspacing)*dY;
imagesc(bws*1e3,deps*1e3,dbah(squeeze(bm_tgc)),[-60 0])
ax=gca;
ax.FontSize=15;
ax.FontWeight='Bold';
colormap gray, colorbar
xlabel('Lateral position [mm]','interpreter','None','fontsize',15,'FontWeight','Bold'),
ylabel('Depth [mm]','interpreter','None','fontsize',15,'FontWeight','Bold')
axis equal, axis tight
%title('bell shaped lesion filled with air bubbles in plain water 5 MHz, 12 ppw, TGC applied', 'Interpreter', 'None', 'FontSize', 10, 'FontWeight', 'normal')
saveas(gcf, 'syn_ntdX_128TGC_bmode.svg', 'svg');
saveas(gcf, 'syn_ntdX_128TGC_bmode.png', 'png');

clear all
disp('beamformer is done')