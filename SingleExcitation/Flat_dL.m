clearvars

lambda0=1;

dLmin=-0.25;
dLmax=0.25;
NdL=201;
dLvec=linspace(dLmin,dLmax,NdL);

Lav=15;

Lopt=0;
        
nperp=12;
dperp=0.8;
gamma=1e1;
    
ratiovec=ones(1,length(dLvec));
Shiftvec=ones(1,length(dLvec));
ueigvat=ones(2*nperp^2,length(dLvec));
ueigvat2=ones(2*nperp^2,length(dLvec));

for iterL=1:NdL

    xvec=setXpos(nperp,dperp);
    yvec=setYpos(nperp,dperp);
    zvec1=-(Lav+dLvec(iterL))*ones(1,nperp^2);
    zvec=[zvec1,-zvec1];

    gvec=ones(size(xvec));
    gvec=gvec/norm(gvec);
    [gammad,gammab,Deltalight,overlaps,Deltabright,ueig,ueig2]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

    ratiovec(iterL)=gammad/gammab;
    Shiftvec(iterL)=(Deltalight-Deltabright)/gamma;
    
    ueigvat(:,iterL)=ueig;
    ueigvat2(:,iterL)=ueig2;
end
%%
figure
plot(Lav+dLvec,ratiovec)