function [] = launch_fullwave2_try6_nln_relaxing4_pzero(c0,omega0,wX,wY,duration,p0,ppw,cfl,cmap,rhomap,Amap,betamap,incoords,incoordszero,outcoords,icmat,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: FEB 5, 2019
% LAST MODIFIED: APR 20, 2021
% Write simulation files
% Add pzero support
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optargin = size(varargin,2);


%nbdy = 20+8; % number of boundary points for PML + stencil
M=8; nbdy=40; 
lambda = c0/omega0*2*pi
				%nX = round(wX/lambda)
nX=size(cmap,1); nY=size(cmap,2);
%nXe = round(wX/lambda*ppw)+2*(nbdy+M);  % number of lateral elements
%nYe = round(wY/lambda*ppw)+2*(nbdy+M);  % number of depth elements
nXe = nX+2*(nbdy+M);  % number of lateral elements
nYe = nY+2*(nbdy+M);  % number of depth elements
nT = round(duration*c0/lambda*ppw/cfl);

dX = c0/omega0*2*pi/ppw
dY = c0/omega0*2*pi/ppw
dT = dY/c0*cfl


%Nmap=(1+boveramap/2)./(rhomap.*cmap.^4);

c=extendMap(cmap,nbdy+M);
rho=extendMap(rhomap,nbdy+M);
beta=extendMap(betamap,nbdy+M);
A=extendMap(Amap,nbdy+M);
K=c.^2.*rho;
size(c);

ncoords = size(incoords,1);
ncoordszero = size(incoordszero,1);
ncoordsout = size(outcoords,1);
incoords = incoords+nbdy+M;
incoordszero = incoordszero+nbdy+M;
outcoords = outcoords+nbdy+M;

omega0_alpha=1e6*2*pi;

nTic=size(icmat,2);
if(optargin)
    nTic=varargin{1};
end



nu=1/2;

r=cfl;



%For 2D modeling:
d(2,1) = -8.74634088067635e-4 * r^7-1.80530560296097e-3 * r^6-4.40512972481673e-4 * r^5 + 4.74018847663366e-3 * r^4-1.93097802254349e-5 * r^3-2.92328221171893e-1 * r^2-6.58101498708345e-8 * r + 1.25420636437969;
d(3,1) = 7.93317828964018e-4 * r^7 + 1.61433256585486e-3 * r^6 + 3.97244786277123e-4 * r^5 + 5.46057645976549e-3 * r^4 + 1.73781972873916e-5 * r^3 + 5.88754971188371e-2 * r^2 + 5.91706982879834e-8 * r-1.23406473759703e-1;
d(4,1) = -6.50217700538851e-4 * r^7-1.16449260340413e-3 * r^6-3.24403734066325e-4 * r^5-9.11483710059994e-3 * r^4-1.41739982312600e-5 * r^3 + 2.33184077551615e-2 * r^2-4.82326094707544e-8 * r + 3.46342451534453e-2;
d(5,1) = 4.67529510541428e-4 * r^7 + 7.32736676632388e-4 * r^6 + 2.32444388955328e-4 * r^5 + 8.46419766685254e-3 * r^4 + 1.01438593426278e-5 * r^3-3.17586249260511e-2 * r^2 + 3.44988852042879e-8 * r-1.19674942518101e-2;
 d(6,1) = -2.98416281187033e-4 * r^7-3.99380750669364e-4 * r^6-1.48203388388213e-4 * r^5-6.01788793192501e-3 * r^4-6.46543538517443e-6 * r^3 + 2.41912754935119e-2 * r^2-2.19855171569984e-8 * r + 4.15554391204146e-3;
d(7,1) = 1.67882669698981e-4 * r^7 + 1.88195874702691e-4 * r^6 + 8.30579218603960e-5 * r^5 + 3.48461963201376e-3 * r^4 + 3.61873162287129e-6 * r^3-1.49875789940005e-2 * r^2 + 1.22979142197165e-8 * r-1.29213888778954e-3;
d(8,1) = -6.22209937489143e-5 * r^7-6.44890425871692e-5 * r^6-3.02936928954918e-5 * r^5-1.33386143898282e-3 * r^4-1.31215186728213e-6 * r^3 + 6.70228205200379e-3 * r^2-4.44653967516776e-9 * r + 3.15659916047599e-4;
d(9,1) = 6.84740881090240e-6 * r^7 + 1.14082245705934e-5 * r^6 + 3.03727593705750e-6 * r^5 + 2.36122782444105e-4 * r^4 + 1.26768491232397e-7 * r^3-1.53347270556276e-3 * r^2 + 4.21617557752767e-10 * r-4.51948990428065e-5;
d(2,2) = 2.13188763071246e-6 * r^7-7.41025068776257e-5 * r^6 + 2.31652037371554e-6 * r^5-2.59495924602038e-3 * r^4 + 1.20637183170338e-7 * r^3 + 5.21123771632193e-2 * r^2 + 4.42258843694177e-10 * r-4.20967682664542e-7;



dmap=zeros(9,2,round(max(max(c)))-round(min(min(c))));
for i=1:round(max(max(c))-min(min(c)))+1
    r=((i-1)+min(min(c)))*dT/dX;
    
    dmap(2,1,i) = -8.74634088067635e-4 * r^7-1.80530560296097e-3 * r^6-4.40512972481673e-4 * r^5 + 4.74018847663366e-3 * r^4-1.93097802254349e-5 * r^3-2.92328221171893e-1 * r^2-6.58101498708345e-8 * r + 1.25420636437969;
dmap(3,1,i) = 7.93317828964018e-4 * r^7 + 1.61433256585486e-3 * r^6 + 3.97244786277123e-4 * r^5 + 5.46057645976549e-3 * r^4 + 1.73781972873916e-5 * r^3 + 5.88754971188371e-2 * r^2 + 5.91706982879834e-8 * r-1.23406473759703e-1;
dmap(4,1,i) = -6.50217700538851e-4 * r^7-1.16449260340413e-3 * r^6-3.24403734066325e-4 * r^5-9.11483710059994e-3 * r^4-1.41739982312600e-5 * r^3 + 2.33184077551615e-2 * r^2-4.82326094707544e-8 * r + 3.46342451534453e-2;
dmap(5,1,i) = 4.67529510541428e-4 * r^7 + 7.32736676632388e-4 * r^6 + 2.32444388955328e-4 * r^5 + 8.46419766685254e-3 * r^4 + 1.01438593426278e-5 * r^3-3.17586249260511e-2 * r^2 + 3.44988852042879e-8 * r-1.19674942518101e-2;
 dmap(6,1,i) = -2.98416281187033e-4 * r^7-3.99380750669364e-4 * r^6-1.48203388388213e-4 * r^5-6.01788793192501e-3 * r^4-6.46543538517443e-6 * r^3 + 2.41912754935119e-2 * r^2-2.19855171569984e-8 * r + 4.15554391204146e-3;
dmap(7,1,i) = 1.67882669698981e-4 * r^7 + 1.88195874702691e-4 * r^6 + 8.30579218603960e-5 * r^5 + 3.48461963201376e-3 * r^4 + 3.61873162287129e-6 * r^3-1.49875789940005e-2 * r^2 + 1.22979142197165e-8 * r-1.29213888778954e-3;
dmap(8,1,i) = -6.22209937489143e-5 * r^7-6.44890425871692e-5 * r^6-3.02936928954918e-5 * r^5-1.33386143898282e-3 * r^4-1.31215186728213e-6 * r^3 + 6.70228205200379e-3 * r^2-4.44653967516776e-9 * r + 3.15659916047599e-4;
dmap(9,1,i) = 6.84740881090240e-6 * r^7 + 1.14082245705934e-5 * r^6 + 3.03727593705750e-6 * r^5 + 2.36122782444105e-4 * r^4 + 1.26768491232397e-7 * r^3-1.53347270556276e-3 * r^2 + 4.21617557752767e-10 * r-4.51948990428065e-5;
dmap(2,2,i) = 2.13188763071246e-6 * r^7-7.41025068776257e-5 * r^6 + 2.31652037371554e-6 * r^5-2.59495924602038e-3 * r^4 + 1.20637183170338e-7 * r^3 + 5.21123771632193e-2 * r^2 + 4.42258843694177e-10 * r-4.20967682664542e-7;
end
r=cfl;
   
dcmap=round(c)-min(min(c))+1;
%% PMLS %%
kpml=1;
dpmlxOld=zeros(1,nXe); alphapmlxOld=dpmlxOld; apmlxOld=dpmlxOld; bplm=dpmlxOld; 
dpmlxOld=zeros(1,nXe); alphapmlxOld=dpmlxOld; apmlxOld=dpmlxOld; bpmlxOld=dpmlxOld+1;
dpmlxOld2=zeros(1,nXe); alphapmlxOld2=dpmlxOld; apmlxOld2=dpmlxOld; bpmlxOld2=dpmlxOld+1;
L=dX*nbdy;
Rc=1e-30;
d0=-3*c0*log(Rc)/(2*L)
for i=1:nbdy
    dpmlxOld(i+(nXe-M-nbdy+1))=d0*(i/nbdy)^2;
    dpmlxOld(M+nbdy+1-i)=d0*(i/nbdy)^2;
end
for i=1:nbdy
    alphapmlxOld(i+(nXe-M-nbdy+1))=omega0/2*(nbdy-i)/nbdy;
    alphapmlxOld(M+nbdy+1-i)=omega0/2*(nbdy-i)/nbdy;
end
alphapmlxOld=alphapmlxOld/10;
%bpmlxOld=exp(-(dpmlxOld/kpml+alphapmlxOld)*dT);
%apmlxOld=dpmlxOld/(kpml*(dpmlxOld+kpml*alphapmlxOld))*(bpmlxOld-1);
[apmlxOld bpmlxOld]=ab(dpmlxOld,kpml,alphapmlxOld,dT);

for i=1:nbdy
    dpmlxOld2(i+(nXe-M-nbdy+1))=d0*((i-1/2)/nbdy)^2;
    dpmlxOld2(M+nbdy+1-i-1)=d0*((i-1/2)/nbdy)^2;
end
for i=1:nbdy
    alphapmlxOld2(i+(nXe-M-nbdy+1))=omega0/2*(nbdy-(i-1/2))/nbdy;
    alphapmlxOld2(M+nbdy+1-i-1)=omega0/2*(nbdy-(i-1/2))/nbdy;
end
alphapmlxOld2=alphapmlxOld2/10;
%bpmlxOld2=exp(-(dpmlxOld2/kpml+alphapmlxOld2)*dT);
%apmlxOld2=dpmlxOld2/(kpml*(dpmlxOld2+kpml*alphapmlxOld2))*(bpmlxOld2-1);
[apmlxOld2 bpmlxOld2]=ab(dpmlxOld2,kpml,alphapmlxOld2,dT);

plot(apmlxOld), hold on
plot(apmlxOld2,'r')
%plot(alphapmlxOld2,'g')
hold off

plot(bpmlxOld), hold on
plot(bpmlxOld2,'r')
hold off


%apmlx=apml; bpmlx=bpml; apmlx2=apml2; bpmlx2=bpml2;
%apmlx(nXe-(M+nbdy)+1:nXe)=fliplr(apmlx(1:M+nbdy));
%bpmlx(nXe-(M+nbdy)+1:nXe)=fliplr(bpmlx(1:M+nbdy));
%apmlx2(nXe-(M+nbdy)+1:nXe)=fliplr(apmlx2(1:M+nbdy));
%bpmlx2(nXe-(M+nbdy)+1:nXe)=fliplr(bpmlx2(1:M+nbdy));

dpmlyOld=zeros(1,nYe); alphapmlyOld=dpmlyOld; apmlyOld=dpmlyOld; bplm=dpmlyOld; 
dpmlyOld=zeros(1,nYe); alphapmlyOld=dpmlyOld; apmlyOld=dpmlyOld; bpmlyOld=dpmlyOld+1; 
dpmlyOld2=zeros(1,nYe); alphapmlyOld2=dpmlyOld; apmlyOld2=dpmlyOld; bpmlyOld2=dpmlyOld+1;

for i=1:nbdy
    dpmlyOld(i+(nYe-M-nbdy+1))=d0*(i/nbdy)^2;
    dpmlyOld(M+nbdy+1-i)=d0*(i/nbdy)^2;
end
for i=1:nbdy
    alphapmlyOld(i+(nYe-M-nbdy+1))=omega0/2*(nbdy-i)/nbdy;
    alphapmlyOld(M+nbdy+1-i)=omega0/2*(nbdy-i)/nbdy;
end
alphapmlyOld=alphapmlyOld/10;
%bpmlyOld=exp(-(dpmlyOld/kpml+alphapmlyOld)*dT);
%apmlyOld=dpmlyOld/(kpml*(dpmlyOld+kpml*alphapmlyOld))*(bpmlyOld-1);
[apmlyOld bpmlyOld]=ab(dpmlyOld,kpml,alphapmlyOld,dT);

for i=1:nbdy
    dpmlyOld2(i+(nYe-M-nbdy+1))=d0*((i-1/2)/nbdy)^2;
    dpmlyOld2(M+nbdy+1-i-1)=d0*((i-1/2)/nbdy)^2;
end
for i=1:nbdy
    alphapmlyOld2(i+(nYe-M-nbdy+1))=omega0/2*(nbdy-(i-1/2))/nbdy;
    alphapmlyOld2(M+nbdy+1-i-1)=omega0/2*(nbdy-(i-1/2))/nbdy;
end
alphapmlyOld2=alphapmlyOld2/10;
%bpmlyOld2=exp(-(dpmlyOld2/kpml+alphapmlyOld2)*dT);
%apmlyOld2=dpmlyOld2/(kpml*(dpmlyOld2+kpml*alphapmlyOld2))*(bpmlyOld2-1);
[apmlyOld2 bpmlyOld2]=ab(dpmlyOld2,kpml,alphapmlyOld2,dT);

%apmlyOld=apml; bpmlyOld=bpml; apmlyOld2=apml2; bpmlyOld2=bpml2;


plot(apmlyOld), hold on
plot(apmlyOld2,'r')
hold off

plot(bpmlyOld), hold on
plot(bpmlyOld2,'r')
hold off

%plot(apmlyOld(M+1:M+nbdy)-apmlyOld(nXe-M+1:-1:nXe-M-nbdy+2))
%plot(apmlyOld2(M:M+nbdy-1)-apmlyOld2(nXe-M+1:-1:nXe-M-nbdy+2))

%pp=zeros(nXe,nYe,nT);


psiu1=zeros(2,nXe,nYe); 
psiu2=zeros(2,nXe,nYe); 
psiw1=zeros(2,nXe,nYe); 
psiw2=zeros(2,nXe,nYe); 
psix1=zeros(2,nXe,nYe); 
psix2=zeros(2,nXe,nYe); 
psiy1=zeros(2,nXe,nYe); 
psiy2=zeros(2,nXe,nYe); 

load xstart_attenuation_dispersion_curves2

nXe
nYe
size(A)

kappax=zeros(nXe,nYe)+1-(1-xstart(1))*sqrt(A/0.5); %ones(nXe,nYe)-0.012;
kappay=zeros(nXe,nYe)+1-(1-xstart(1))*sqrt(A/0.5); %ones(nXe,nYe)-0.012;
kappau=zeros(nXe,nYe)+1-(1-xstart(6))*sqrt(A/0.5); %ones(nXe,nYe)-0.015;
kappaw=zeros(nXe,nYe)+1-(1-xstart(6))*sqrt(A/0.5); %ones(nXe,nYe)-0.015;

[apmlx1 bpmlx1]=ab(xstart(2),xstart(1),xstart(3),dT./A*0.5);%1.0278e-11);
[apmlx2 bpmlx2]=ab(xstart(4),xstart(1),xstart(5),dT./A*0.5);
[apmlu1 bpmlu1]=ab(xstart(7),xstart(6),xstart(8),dT./A*0.5);
[apmlu2 bpmlu2]=ab(xstart(9),xstart(6),xstart(10),dT./A*0.5);



if(0)
    
    [apmlx1 bpmlx1]=ab(xstart(2),xstart(1),xstart(3),dT*omega0/omega0_alpha./A*0.5);%1.0278e-11);
[apmlx2 bpmlx2]=ab(xstart(4),xstart(1),xstart(5),dT*omega0/omega0_alpha./A*0.5);
[apmlu1 bpmlu1]=ab(xstart(7),xstart(6),xstart(8),dT*omega0/omega0_alpha./A*0.5);
[apmlu2 bpmlu2]=ab(xstart(9),xstart(6),xstart(10),dT*omega0/omega0_alpha./A*0.5);


[apmlx1 bpmlx1]=ab(xstart(2),xstart(1),xstart(3),dT/3.8894e+06);%1.0278e-11);
[apmlx2 bpmlx2]=ab(xstart(4),xstart(1),xstart(5),dT/3.8894e+06);
[apmlu1 bpmlu1]=ab(xstart(7),xstart(6),xstart(8),dT/3.8894e+06);
[apmlu2 bpmlu2]=ab(xstart(9),xstart(6),xstart(10),dT/3.8894e+06);
end
apmly1=apmlx1;
bpmly1=bpmlx1;
apmly2=apmlx2;
bpmly2=bpmlx2;
apmlw1=apmlu1;
bpmlw1=bpmlu1;
apmlw2=apmlu2;
bpmlw2=bpmlu2;

apmlu1=zeros(nXe,nYe)+apmlu1;
bpmlu1=zeros(nXe,nYe)+bpmlu1;
apmlw1=zeros(nXe,nYe)+apmlw1;
bpmlw1=zeros(nXe,nYe)+bpmlw1;
apmlx1=zeros(nXe,nYe)+apmlx1;
bpmlx1=zeros(nXe,nYe)+bpmlx1;
apmly1=zeros(nXe,nYe)+apmly1;
bpmly1=zeros(nXe,nYe)+bpmly1;

apmlu2=zeros(nXe,nYe)+apmlu2;
bpmlu2=zeros(nXe,nYe)+bpmlu2;
apmlw2=zeros(nXe,nYe)+apmlw2;
bpmlw2=zeros(nXe,nYe)+bpmlw2;
apmlx2=zeros(nXe,nYe)+apmlx2;
bpmlx2=zeros(nXe,nYe)+bpmlx2;
apmly2=zeros(nXe,nYe)+apmly2;
bpmly2=zeros(nXe,nYe)+bpmly2;


% localize PML region
pmlmaskx=zeros(nXe,nYe);
pmlmasky=zeros(nXe,nYe);
for i=1:nbdy
    pmlmaskx(i+(nXe-M-nbdy+1),:)=(i/nbdy);
    pmlmaskx(M+nbdy+1-i,:)=(i/nbdy);
    pmlmasky(:,i+(nYe-M-nbdy+1))=(i/nbdy);
    pmlmasky(:,M+nbdy+1-i)=(i/nbdy);
end
pmlmaskx(1:M,:)=1; pmlmaskx(nXe-M+1:nXe,:)=1;
pmlmasky(:,1:M)=1; pmlmasky(:,nYe-M+1:nYe,:)=1;



plot(apmlxOld), hold on
plot(apmlxOld2,'r')
plot(apmlx1(:,round(nYe/2)),'g')
hold off

plot(bpmlxOld), hold on
plot(bpmlxOld2,'r')
plot(bpmlx1(:,round(nYe/2)),'g')
hold off
 
imagesc(apmlx1.*(1-pmlmaskx)+(apmlx1+apmlu1)./2.*pmlmaskx)

[apmlx1 apmlu1]=pml_gradient_mask2(nXe,nYe,nbdy/2,M+nbdy/2,1,apmlx1,apmlu1);
[apmlx2 apmlu2]=pml_gradient_mask2(nXe,nYe,nbdy/2,M+nbdy/2,1,apmlx2,apmlu2);
[apmly1 apmlw1]=pml_gradient_mask2(nXe,nYe,nbdy/2,M+nbdy/2,2,apmly1,apmlw1);
[apmly2 apmlw2]=pml_gradient_mask2(nXe,nYe,nbdy/2,M+nbdy/2,2,apmly2,apmlw2);
[bpmlx1 bpmlu1]=pml_gradient_mask2(nXe,nYe,nbdy/2,M+nbdy/2,1,bpmlx1,bpmlu1);
[bpmlx2 bpmlu2]=pml_gradient_mask2(nXe,nYe,nbdy/2,M+nbdy/2,1,bpmlx2,bpmlu2);
[bpmly1 bpmlw1]=pml_gradient_mask2(nXe,nYe,nbdy/2,M+nbdy/2,2,bpmly1,bpmlw1);
[bpmly2 bpmlw2]=pml_gradient_mask2(nXe,nYe,nbdy/2,M+nbdy/2,2,bpmly2,bpmlw2);

[apmlx1 tmp]=pml_gradient_mask2(nXe,nYe,nbdy,M,1,apmlx1,apmlxOld'*ones(1,nYe));
[apmlu1 tmp]=pml_gradient_mask2(nXe,nYe,nbdy,M,1,apmlu1,apmlxOld2'*ones(1,nYe));
[apmly1 tmp]=pml_gradient_mask2(nXe,nYe,nbdy,M,2,apmly1,ones(nXe,1)*apmlyOld);
[apmlw1 tmp]=pml_gradient_mask2(nXe,nYe,nbdy,M,2,apmlw1,ones(nXe,1)*apmlyOld2);
[bpmlx1 tmp]=pml_gradient_mask2(nXe,nYe,nbdy,M,1,bpmlx1,bpmlxOld'*ones(1,nYe));
[bpmlu1 tmp]=pml_gradient_mask2(nXe,nYe,nbdy,M,1,bpmlu1,bpmlxOld2'*ones(1,nYe));
[bpmly1 tmp]=pml_gradient_mask2(nXe,nYe,nbdy,M,2,bpmly1,ones(nXe,1)*bpmlyOld);
[bpmlw1 tmp]=pml_gradient_mask2(nXe,nYe,nbdy,M,2,bpmlw1,ones(nXe,1)*bpmlyOld2);

if(0)
    
[apmlx2 tmp]=pml_gradient_mask(nXe,nYe,nbdy,M,1,apmlx2,apmlxOld'*ones(1,nYe));
[apmlu2 tmp]=pml_gradient_mask(nXe,nYe,nbdy,M,1,apmlu2,apmlxOld2'*ones(1,nYe));
[apmly2 tmp]=pml_gradient_mask(nXe,nYe,nbdy,M,2,apmly2,ones(nXe,1)*apmlyOld);
[apmlw2 tmp]=pml_gradient_mask(nXe,nYe,nbdy,M,2,apmlw2,ones(nXe,1)*apmlyOld2);
[bpmlx2 tmp]=pml_gradient_mask(nXe,nYe,nbdy,M,1,bpmlx2,bpmlxOld'*ones(1,nYe));
[bpmlu2 tmp]=pml_gradient_mask(nXe,nYe,nbdy,M,1,bpmlu2,bpmlxOld2'*ones(1,nYe));
[bpmly2 tmp]=pml_gradient_mask(nXe,nYe,nbdy,M,2,bpmly2,ones(nXe,1)*bpmlyOld);
[bpmlw2 tmp]=pml_gradient_mask(nXe,nYe,nbdy,M,2,bpmlw2,ones(nXe,1)*bpmlyOld2);


plot(bpmlx1(:,round(nYe/2)),'c'), hold on
plot(bpmlxOld)
plot(apmlx1(:,round(nYe/2)),'g'), hold on
plot(apmlxOld,'r')
hold off

plot(bpmlx2(:,round(nYe/2)),'c'), hold on
plot(bpmlxOld2)
plot(apmlx2(:,round(nYe/2)),'g'), hold on
plot(apmlxOld2,'r')
hold off
end


fid = fopen(['c.dat'],'wb'); fwrite(fid,c','float'); fclose(fid);
fid = fopen(['K.dat'],'wb'); fwrite(fid,K','float'); fclose(fid);
fid = fopen(['rho.dat'],'wb'); fwrite(fid,rho','float'); fclose(fid);
fid = fopen(['beta.dat'],'wb'); fwrite(fid,beta','float'); fclose(fid);

writeCoords('icc.dat',incoords-1);
writeCoords('icczero.dat',incoordszero-1);
writeCoords('outc.dat',outcoords-1);
writeIC('icmat.dat',icmat');

writeVabs('float',dX,'dX',dY,'dY',dT,'dT',c0,'c0');
writeVabs('int',nXe,'nX',nYe,'nY',nT,'nT',ncoords,'ncoords',ncoordszero,'ncoordszero',ncoordsout,'ncoordsout',nTic,'nTic');


fid = fopen(['d.dat'],'wb'); fwrite(fid,d','float'); fclose(fid);

fid = fopen(['dmap.dat'],'wb'); fwrite(fid,zeros(0),'float'); fclose(fid);

fid = fopen(['dmap.dat'],'ab'); 
fwrite(fid,zeros(0),'float'); 
for i=1:size(dmap,1)
    for j=1:size(dmap,2)
        for k=1:size(dmap,3)
            fwrite(fid,dmap(i,j,k),'float');
        end
    end
end
fclose(fid);

ndmap=size(dmap,3)
writeVabs('int',ndmap,'ndmap');

fid = fopen(['dcmap.dat'],'wb'); fwrite(fid,(dcmap-1)','int'); fclose(fid);

if(0)
    fid = fopen(['dpmlx'],'wb'); fwrite(fid,dpmlx,'float'); fclose(fid);
    fid = fopen(['alphapmlx'],'wb'); fwrite(fid,alphapmlx,'float'); fclose(fid);
    fid = fopen(['bpmlx'],'wb'); fwrite(fid,bpmlx,'float'); fclose(fid);
    fid = fopen(['apmlx'],'wb'); fwrite(fid,apmlx,'float'); fclose(fid);
    fid = fopen(['dpmlx2'],'wb'); fwrite(fid,dpmlx2,'float'); fclose(fid);
    fid = fopen(['alphapmlx2'],'wb'); fwrite(fid,alphapmlx2,'float'); fclose(fid);
    fid = fopen(['bpmlx2'],'wb'); fwrite(fid,bpmlx2,'float'); fclose(fid);
    fid = fopen(['apmlx2'],'wb'); fwrite(fid,apmlx2,'float'); fclose(fid);
    fid = fopen(['dpmly'],'wb'); fwrite(fid,dpmly,'float'); fclose(fid);
    fid = fopen(['alphapmly'],'wb'); fwrite(fid,alphapmly,'float'); fclose(fid);
    fid = fopen(['bpmly'],'wb'); fwrite(fid,bpmly,'float'); fclose(fid);
    fid = fopen(['apmly'],'wb'); fwrite(fid,apmly,'float'); fclose(fid);
    fid = fopen(['dpmly2'],'wb'); fwrite(fid,dpmly2,'float'); fclose(fid);
    fid = fopen(['alphapmly2'],'wb'); fwrite(fid,alphapmly2,'float'); fclose(fid);
    fid = fopen(['bpmly2'],'wb'); fwrite(fid,bpmly2,'float'); fclose(fid);
    fid = fopen(['apmly2'],'wb'); fwrite(fid,apmly2,'float'); fclose(fid);
    
end
 

  fid = fopen(['kappax.dat'],'wb'); fwrite(fid,kappax','float'); fclose(fid);
  fid = fopen(['kappay.dat'],'wb'); fwrite(fid,kappay','float'); fclose(fid);
  fid = fopen(['kappau.dat'],'wb'); fwrite(fid,kappau','float'); fclose(fid);
  fid = fopen(['kappaw.dat'],'wb'); fwrite(fid,kappaw','float'); fclose(fid);
  
  fid = fopen(['apmlu1.dat'],'wb'); fwrite(fid,apmlu1','float'); fclose(fid);
  fid = fopen(['bpmlu1.dat'],'wb'); fwrite(fid,bpmlu1','float'); fclose(fid);
  fid = fopen(['apmlw1.dat'],'wb'); fwrite(fid,apmlw1','float'); fclose(fid);
  fid = fopen(['bpmlw1.dat'],'wb'); fwrite(fid,bpmlw1','float'); fclose(fid);
  fid = fopen(['apmlx1.dat'],'wb'); fwrite(fid,apmlx1','float'); fclose(fid);
  fid = fopen(['bpmlx1.dat'],'wb'); fwrite(fid,bpmlx1','float'); fclose(fid);
  fid = fopen(['apmly1.dat'],'wb'); fwrite(fid,apmly1','float'); fclose(fid);
  fid = fopen(['bpmly1.dat'],'wb'); fwrite(fid,bpmly1','float'); fclose(fid);
  
  fid = fopen(['apmlu2.dat'],'wb'); fwrite(fid,apmlu2','float'); fclose(fid);
  fid = fopen(['bpmlu2.dat'],'wb'); fwrite(fid,bpmlu2','float'); fclose(fid);
  fid = fopen(['apmlw2.dat'],'wb'); fwrite(fid,apmlw2','float'); fclose(fid);
  fid = fopen(['bpmlw2.dat'],'wb'); fwrite(fid,bpmlw2','float'); fclose(fid);
  fid = fopen(['apmlx2.dat'],'wb'); fwrite(fid,apmlx2','float'); fclose(fid);
  fid = fopen(['bpmlx2.dat'],'wb'); fwrite(fid,bpmlx2','float'); fclose(fid);
  fid = fopen(['apmly2.dat'],'wb'); fwrite(fid,apmly2','float'); fclose(fid);
  fid = fopen(['bpmly2.dat'],'wb'); fwrite(fid,bpmly2','float'); fclose(fid);


