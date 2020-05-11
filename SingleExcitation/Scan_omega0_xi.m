clearvars

tic

lambda0=1;

omega0min=0.5*lambda0;
omega0max=5*lambda0;
Nomega0=100;
omega0vec=linspace(omega0min,omega0max,Nomega0);

ximin=-2;
ximax=0;
Nxi=20;
xivec=10.^linspace(ximin,ximax,Nxi);

ratiovec=zeros(1,Nxi);
omega0opt=zeros(1,Nxi);
xiopt=zeros(1,Nxi);
        
nperp=8;
dperp=0.7;
gamma=1e1;
        
for iterL=1:Nxi

    L=0.5*xivec(iterL)*(nperp*dperp)^2/lambda0;
    
    ratiovec(iterL)=1;
    
    for iteromega0=1:Nomega0
        
        omega0=omega0vec(iteromega0);
        zR=pi*omega0^2/lambda0;

        xvec=setXpos(nperp,dperp);
        yvec=setYpos(nperp,dperp);
        zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

        gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
        [gammad,gammab,Deltalight,overlaps]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

        if gammad/gammab<ratiovec(iterL)
            ratiovec(iterL)=gammad/gammab;
            omega0opt(iterL)=omega0;
            xiopt(iterL)=2*max(zvec)*lambda0/(nperp*dperp)^2;
        end

    end

end
toc

figure
plot(xiopt,omega0opt)

figure
loglog(xiopt,ratiovec)
