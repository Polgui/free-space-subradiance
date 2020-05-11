function omega0 = find_omega0opt(lambda0,nperp,dperp,L)
na=nperp^2;
gamma=1;
    omega0min=0.5*lambda0;
    omega0max=5*lambda0;
    Nomega0=100;
    omega0vec=linspace(omega0min,omega0max,Nomega0);

    ratiovec=zeros(1,Nomega0);

    minzvec=zeros(1,Nomega0);
    maxzvec=zeros(1,Nomega0);
    allzvec=zeros(2*na,Nomega0);

    n0vec=zeros(1,Nomega0);

    for iteromega0=1:Nomega0

        omega0=omega0vec(iteromega0);
        zR=pi*omega0^2/lambda0;


        xvec=setXpos(nperp,dperp);
        yvec=setYpos(nperp,dperp);
        zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

        gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
        [gammad,gammab,Deltalight,overlaps]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

        ratiovec(iteromega0)=gammad/gammab;

        minzvec(iteromega0)=min(abs(zvec));
        maxzvec(iteromega0)=max(abs(zvec));

        allzvec(:,iteromega0)=abs(zvec);


    end
    
[~,iter]=min(ratiovec);

omega0=omega0vec(iter);

end

