clearvars

lambda0=1;

nperp=16;
dperp=0.5;
gamma=1e1;

ximin=-1;
ximax=1;
Nxi=20;
xivec=10.^linspace(ximin,ximax,Nxi);
ximinvec=zeros(1,length(xivec));

ratiominvec=ones(1,length(xivec));

for iterxi=1:Nxi
    iterxi/Nxi

Lav=0.5*xivec(iterxi)*(nperp*dperp)^2/lambda0;

dLmin=-0.125;
dLmax=0.125;
NdL=101;
dLvec=linspace(dLmin,dLmax,NdL);

ratiovec=ones(1,length(dLvec));

for iterL=1:NdL

    xvec=setXpos(nperp,dperp);
    yvec=setYpos(nperp,dperp);
    zvec1=-(Lav+dLvec(iterL))*ones(1,nperp^2);
    zvec=[zvec1,-zvec1];

    gvec=ones(size(xvec));
    gvec=gvec/norm(gvec);
    [gammad,gammab,Deltalight,overlaps,Deltabright,ueig,ueig2]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

    ratiovec(iterL)=gammad/gammab;
end
[ratiominvec(iterxi),ind]=min(ratiovec);
ximinvec(iterxi)=2*(Lav+dLvec(ind))*lambda0/(nperp*dperp)^2;

end

%%

ratiooptvec=ones(1,length(xivec));
xioptvec=zeros(1,length(xivec));

for iterxi=1:Nxi
iterxi/Nxi
L=0.5*xivec(iterxi)*(nperp*dperp)^2/lambda0;

omega0=find_omega0opt(lambda0,nperp,dperp,L);
zR=pi*omega0^2/lambda0;

xvec=setXpos(nperp,dperp);
yvec=setYpos(nperp,dperp);
[zvec,n0]=setZpos(nperp,L,lambda0,zR,xvec,yvec,lambda0/4);

gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
[gammad,gammab,Deltalight,overlaps,Deltabright,ueig,ueig2,lam,u]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

gammad;
gammab;
ratiooptvec(iterxi)=gammad/gammab;
xioptvec(iterxi)=2*(lambda0/4)*n0*lambda0/(nperp*dperp)^2;
end

%%
figure
loglog(ximinvec,ratiominvec)
hold on
loglog(xioptvec,ratiooptvec)

%%

figure
loglog(xivec,0.22*xivec.^(1.25))
