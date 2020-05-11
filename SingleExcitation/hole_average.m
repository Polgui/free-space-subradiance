%% gammad/gammab
clearvars

lambda0=1;

nperp=12;
dperp=0.5;
L=5;
gamma=1e1;

na=nperp^2;

omega0=find_omega0opt(lambda0,nperp,dperp,L);
zR=pi*omega0^2/lambda0;
%%
gammab_vec=zeros(1,2*na);
gammad_vec=zeros(1,2*na);
gvec_mean=0;

for iterhole=1:2*na
    iterhole

    xvec=setXpos(nperp,dperp);
    yvec=setYpos(nperp,dperp);
    zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

    gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
    gvec_mean=gvec_mean+abs(gvec(iterhole))^2;

    xvec(iterhole)=(nperp*dperp)*10000;

    gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
    [gammad,gammab,Deltalight,overlaps,Deltabright,ueig]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

    gammad_vec(iterhole)=gammad;
    gammab_vec(iterhole)=gammab;
end

%% Plots

figure
plot(gammab_vec)
hold on
plot(gammad_vec)
plot(gammad_vec./gammab_vec)

mean(gammad_vec)/mean(gammab_vec)