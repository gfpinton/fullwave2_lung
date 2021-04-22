function [pxducer2] = gen_txrx_pxducer(basedir,ncoordsout,idc,nT,nTxRx,xdcReduction,flt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2020-07-22
% LAST MODIFIED: 2020-07-22
% Generate transducer data from multiple txrx sims
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basedir : directory with all the txrx events
% nT : number of time points
% xdcReduction : label vector to average sim points to element points should have values in between 1 and number of transducer elements
% flt : filter vector of length nT

  if (length(nTxRx)==1)
    nTxRxvec=1:nTxRx;
  else
    nTxRxvec=nTxRx;
  end
  
basedir
pxducer2=zeros(nT,max(xdcReduction),nTxRx,'single');
for n=nTxRxvec
  outdir=[basedir '/txrx_' num2str(n)]; 
  nRun=sizeOfFile([outdir '/genout.dat'])/4/ncoordsout;
  if(nRun>nT)
    disp(['Warning, nRun>nT, nRun=' num2str(nRun) ' nT=' num2str(nT)]);
    nRun=nT;
  end
  if(nRun~=0)
    pxducer=readGenoutSlice([outdir '/genout.dat'],0:nRun-2,ncoordsout,idc);
    pxducer(end+1:nT,:)=0; % pad with zeros if sim is not done running
  end
  if(nRun==0)
    pxducer=zeros(nT,ncoordsout);% pad with zeros if sim is not done running
    disp(['Warning, file in ' num2str(n) ' is empty!'])
  end
  pxducer1=ifft(fft(pxducer).*(flt.*ones(1,size(pxducer,2))),'symmetric'); % filter around transmit freq
  for tt=1:max(xdcReduction)
    pxducer2(:,tt,n)=mean(pxducer1(:,find(xdcReduction==tt)),2);
  end
end
  
   
