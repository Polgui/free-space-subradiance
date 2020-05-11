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

%% any number of holes

Nreal=500;

pdef=1e-2;
pdef=2*0.0035;
pdef=0.5;

gammad=0;
gammad_over=0;
gammab_over=0;
gammab=0;
ndefect=0;
ov=0;

xvec=setXpos(nperp,dperp);
yvec=setYpos(nperp,dperp);
zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
[gammad_free]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

gvecB=gvec;
gvecD=gvec.*(-1).^(zvec>0);
    
for iterreal=1:Nreal

    is_defect=binornd(1,pdef,1,2*na);
    is_there=1-is_defect;
    ndefect=ndefect+sum(is_defect)/Nreal;

    gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
    Hnh=zeros(length(xvec),length(xvec));

    for j=1:length(xvec)
        Hnh(j,j)=Hnh(j,j)-is_there(j)*1j*gamma/2.;
        for k=j+1:length(xvec)
           Jhop=is_there(j)*is_there(k)*1.5*gamma*lambda0*GreensTensor(lambda0,sqrt((xvec(k)-xvec(j))^2+(yvec(k)-yvec(j))^2),zvec(k)-zvec(j));
           Hnh(k,j)=Hnh(k,j)-Jhop;
           Hnh(j,k)=Hnh(j,k)-Jhop;
        end
    end

    [u,v]=eig(Hnh);
    lam=diag(v);

    g1=gvec;
    g1(length(xvec)/2+1:end)=0;

    overlaps=conj(g1)*u;
    [~,index]=sort(abs(overlaps).^2);
    maxoverlaps=[lam(index(end)),lam(index(end-1))];

    ueig=u(:,index(end));

    gammab_temp=max(-2*imag(maxoverlaps));
    [gammad_temp,indexgammad]=min(-2*imag(maxoverlaps));

    overlaps=overlaps(index);
    gammad=gammad+gammad_temp/Nreal;
    gammab=gammab+gammab_temp/Nreal;
    ov=ov+2*sum(abs(overlaps(end-1:end)).^2)/Nreal;
    
    gammad_over=gammad_over-2*imag(gvecD*Hnh*gvecD'/Nreal);
    gammab_over=gammab_over-2*imag(gvecB*Hnh*gvecB'/Nreal);
end
ndefect
gammad_free
gammad
gammad_over
gamma*pdef
gammad_free+pdef*(gamma-2*gammad_free)+pdef^2*(gammad_free-gamma)
(gammad_free+pdef*(gamma-2*gammad_free)+pdef^2*(gammad_free-gamma))/
ov
ratiovec=gammad/gammab

%% force only single hole effects

Nreal=200;

pdef=1e-2;
pdef=0.0035;

gammad=0;
gammab=0;
ndefect=0;
ov=0;

for iterreal=1:Nreal
    
    zR=pi*omega0^2/lambda0;

    is_defect=binornd(1,pdef,1,2*na);
    
    ndef=sum(is_defect);
    
    ndefect=ndefect+sum(is_defect)/Nreal;

    for iterdef=1:length(is_defect)
        if is_defect(iterdef)>0.5
            xvec=setXpos(nperp,dperp);
            yvec=setYpos(nperp,dperp);
            zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);
            xvec(iterdef)=xvec(iterdef)+lambda0*10000;

            gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
            [gammad_temp,gammab_temp,Deltalight,overlaps,Deltabright,ueig]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

            gammad=gammad+gammad_temp/(Nreal);
            gammab=gammab+gammab_temp/(Nreal);
            ov=ov+2*sum(abs(overlaps(end-1:end)).^2)/(ndef*Nreal);
        end
    end
end
ndefect
gammad
gamma*pdef
ov
ratiovec=gammad/gammab