%% gammad/gammab
clearvars

nperpvec=[4,8,12,16,20];
Nnperp=length(nperpvec);

pdefvec=10.^linspace(-1,-2,20);
Ndef=length(pdefvec);

Nreal=100;

mycluster=parcluster();
mypool=parpool(mycluster);  

lambda0=1;

dperp=0.8;
L=10;
gamma=1e1;

gammab_av=zeros(Nnperp,Ndef);
gammad_av=zeros(Nnperp,Ndef);
overlapmin=zeros(Nnperp,Ndef);
overlapmax=zeros(Nnperp,Ndef);

for iternperp=1:Nnperp
    nperp=nperpvec(iternperp);
    na=nperp^2;

    omega0=find_omega0opt(lambda0,nperp,dperp,L);
    zR=pi*omega0^2/lambda0;
    disp(nperp)

    for iterdef=1:Ndef
        
        pdef=pdefvec(iterdef);
        
        gammadvec=zeros(1,Nreal);
        gammabvec=zeros(1,Nreal);
        ominvec=zeros(1,Nreal);
        omaxvec=zeros(1,Nreal);
        
        parfor iterreal=1:Nreal

            is_defect=binornd(1,pdef,1,2*na);
            is_there=1-is_defect;
            ndefect=sum(is_defect);

            xvec=setXpos(nperp,dperp);
            yvec=setYpos(nperp,dperp);
            zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

            gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
            Hnh=zeros(length(xvec),length(xvec));

            for j=1:length(xvec)
                Hnh(j,j)=Hnh(j,j)-is_there(j)*1j*gamma/2.;
                for k=j+1:length(xvec)
                   Jhop=1.5*gamma*is_there(j)*is_there(k)*lambda0*GreensTensor(lambda0,sqrt((xvec(k)-xvec(j))^2+(yvec(k)-yvec(j))^2),zvec(k)-zvec(j));
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

            gammab=max(-2*imag(maxoverlaps));
            [gammad,indexgammad]=min(-2*imag(maxoverlaps));

            Deltalight=real(maxoverlaps(indexgammad));

            Deltabright=real(maxoverlaps(3-indexgammad));

            overlaps=overlaps(index);

            gammadvec(iterreal)=gammad;
            gammabvec(iterreal)=gammab;
            sortoverlaps=sort(abs(overlaps).^2);
            ominvec(iterreal)=2*sum(sortoverlaps(1:end-2));
            omaxvec(iterreal)=2*sum(sortoverlaps(end-1:end));
        end
        gammad_av(iternperp,iterdef)=mean(gammadvec);
        gammab_av(iternperp,iterdef)=mean(gammabvec);
        overlapmin(iternperp,iterdef)=mean(ominvec);
        overlapmax(iternperp,iterdef)=mean(omaxvec);
    end
end

delete(gcp)

save('defects_scan.mat','gammad_av','gammab_av','overlapmin','overlapmax')

%%

nperpvec=[4,8,12,16,20];
Nnperp=length(nperpvec);

gammabvec=zeros(1,Nnperp);
gammadvec=zeros(1,Nnperp);

for iternperp=1:Nnperp
    
    dperp=0.5;
    gamma=10;
    lambda0=1;
    nperp=nperpvec(iternperp)
    na=nperp^2;
    L=15;
    omega0=find_omega0opt(lambda0,nperp,dperp,L);

    zR=pi*omega0^2/lambda0;

    xvec=setXpos(nperp,dperp);
    yvec=setYpos(nperp,dperp);
    zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

    gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);

    [gammadvec(iternperp),gammabvec(iternperp)]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);
 
end
%%
save('defects_gvec_mean.mat','gammabvec','gammadvec')

%%

load('defects_scan.mat')
load('defects_gvec_mean.mat')

myBlue=[0.6,0.8,1];

nperpvec=[4,8,12,16,20];
Nnperp=length(nperpvec);

pdefvec=10.^linspace(-3,-1,20);
Ndef=length(pdefvec);

lambda0=1;

dperp=0.5;
L=15;
gamma=1e1;

figure
for iternperp=1:Nnperp
    h=loglog(pdefvec,gammad_av(iternperp,:)./gammab_av(iternperp,:));
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    h.LineWidth=2.0;
    h.MarkerSize=25;
    hold on
    
    gammad=gammadvec(iternperp);
    gammab=gammabvec(iternperp);
    gammad_fit=gammad*(1-2*pdefvec)+pdefvec*gamma+pdefvec.^2*(gammad-gamma);
    gammab_fit=gammab*(1-2*pdefvec)+pdefvec*gamma+pdefvec.^2*(gammad-gamma);;
    
    h=loglog(pdefvec,gammad_fit./gammab_fit,'--');
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    h.LineWidth=2.0;
    h.MarkerSize=25;
    hold on
end

set(gca,'FontSize',22)
ax = gca;
ax.FontName = 'LaTeX';
grid on
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.ZLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
xlabel('$p_d$')
ylabel('$\gamma_d/\gamma_b$')
ax.XLim=[1e-3,1e-1];

%%

load('defects_scan.mat')
load('defects_gvec_mean.mat')

myBlue=[0.6,0.8,1];

nperpvec=[4,8,12,16,20];
Nnperp=length(nperpvec);

pdefvec=10.^linspace(-3,-1,20);
Ndef=length(pdefvec);

lambda0=1;

dperp=0.5;
L=15;
gamma=1e1;

figure
for iternperp=1:Nnperp
    h=plot(pdefvec,exp(-sqrt(gammad_av(iternperp,:)./gammab_av(iternperp,:))));
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    h.LineWidth=2.0;
    h.MarkerSize=25;
    hold on
    
    gammad=gammadvec(iternperp);
    gammab=gammabvec(iternperp);
    gammad_fit=gammad*(1-2*pdefvec)+2*pdefvec*gamma;
    gammab_fit=gammab*(1-2*pdefvec)+2*pdefvec*gamma;
    
    h=plot(pdefvec,exp(-sqrt(gammad_fit./gammab_fit)),'--');
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    h.LineWidth=2.0;
    h.MarkerSize=25;
    hold on
end
