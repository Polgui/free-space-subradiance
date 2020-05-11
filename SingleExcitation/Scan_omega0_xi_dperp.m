clearvars

tic

lambda0=1;

omega0min=0.5*lambda0;
omega0max=5*lambda0;
Nomega0=100;
omega0vec=linspace(omega0min,omega0max,Nomega0);

ximin=-1;
ximax=0;
Nxi=10;
xivec=10.^linspace(ximin,ximax,Nxi);

dperpmin=0.1;
dperpmax=1.0;
Ndperp=10;
dperpvec=linspace(dperpmin,dperpmax,Ndperp);

ratiovec=zeros(Ndperp,Nxi);
omega0opt=zeros(Ndperp,Nxi);
xiopt=zeros(Ndperp,Nxi);
        
nperp=8;
gamma=1e1;

for iterdperp=1:Ndperp
    
    dperp=dperpvec(iterdperp)
        
    for iterL=1:Nxi

        L=0.5*xivec(iterL)*(nperp*dperp)^2/lambda0;

        ratiovec(iterdperp,iterL)=1;

        for iteromega0=1:Nomega0

            omega0=omega0vec(iteromega0);
            zR=pi*omega0^2/lambda0;

            xvec=setXpos(nperp,dperp);
            yvec=setYpos(nperp,dperp);
            zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

            gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
            [gammad,gammab,Deltalight,overlaps]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

            if gammad/gammab<ratiovec(iterdperp,iterL)
                ratiovec(iterdperp,iterL)=gammad/gammab;
                omega0opt(iterdperp,iterL)=omega0;
                xiopt(iterdperp,iterL)=2*max(zvec)*lambda0/(nperp*dperp)^2;
            end

        end

    end
end
toc

figure
for iterdperp=1:Ndperp
    plot(xiopt(iterdperp,:),omega0opt(iterdperp,:))
    hold on
end

figure
for iterdperp=1:Ndperp
    loglog(xiopt(iterdperp,:),ratiovec(iterdperp,:))
    hold on
end
