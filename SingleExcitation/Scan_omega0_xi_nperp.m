clearvars

dperp=.5;

lambda0=1;

omega0min=0.5*lambda0;
omega0max=5*lambda0;
Nomega0=100;
omega0vec=linspace(omega0min,omega0max,Nomega0);

ximin=-1;
ximax=0;
Nxi=20;
xivec=10.^linspace(ximin,ximax,Nxi);

nperpvec=[4,8,12,16,20];
Nnperp=length(nperpvec);

ratiovec=zeros(Nnperp,Nxi);
omega0opt=zeros(Nnperp,Nxi);
xiopt=zeros(Nnperp,Nxi);
        
gamma=1e1;

for iternperp=1:Nnperp
    
    nperp=nperpvec(iternperp)
        
    for iterL=1:Nxi

        L=0.5*xivec(iterL)*(nperp*dperp)^2/lambda0;

        ratiovec(iternperp,iterL)=1;

        for iteromega0=1:Nomega0

            omega0=omega0vec(iteromega0);
            zR=pi*omega0^2/lambda0;

            xvec=setXpos(nperp,dperp);
            yvec=setYpos(nperp,dperp);
            zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

            gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
            [gammad,gammab,Deltalight,overlaps]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

            if gammad/gammab<ratiovec(iternperp,iterL)
                ratiovec(iternperp,iterL)=gammad/gammab;
                omega0opt(iternperp,iterL)=omega0;
                xiopt(iternperp,iterL)=2*max(zvec)*lambda0/(nperp*dperp)^2;
            end

        end

    end
end

figure
for iternperp=1:Nnperp
    plot(xiopt(iternperp,:),omega0opt(iternperp,:))
    hold on
end

figure
set(gcf,'colormap',hot(32))
for iternperp=1:Nnperp
    loglog(xiopt(iternperp,:),ratiovec(iternperp,:))
    hold on
end
