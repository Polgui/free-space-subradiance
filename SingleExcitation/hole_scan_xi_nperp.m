%% gammad/gammab
clearvars

load('Compress.mat')

lambda0=1;

xiopt=xiopt_cell{5};
omega0opt=omega0opt_cell{5};

Nxi=size(xiopt,2);
dperp=0.5;

gamma=1e1;

nperpvec=[4,8,12,16,20];
Nnperp=length(nperpvec);

gammab_av=zeros(Nnperp,Nxi);
gammad_av=zeros(Nnperp,Nxi);

mycluster=parcluster();
mypool=parpool(mycluster);  

for iternperp=1:Nnperp
    
    nperp=nperpvec(iternperp);
    na=nperp^2;
    
    for iterL=1:Nxi

        if 0.5*xiopt(iternperp,iterL)*(nperp*dperp)^2 > 1

            L=0.5*xiopt(iternperp,iterL)*(nperp*dperp)^2/lambda0+0.1;

            omega0=omega0opt(iternperp,iterL);

            zR=pi*omega0^2/lambda0;
            
            gammab_vec=zeros(1,na);
            gammad_vec=zeros(1,na);

            xvec=setXpos(nperp,dperp);
            yvec=setYpos(nperp,dperp);
            zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

            gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
            
            
            parfor iterhole=na+1:2*na

                xvec_hole=xvec;
                xvec_hole(iterhole)=(nperp*dperp)*10000;

                [gammad,gammab]=getShift(gvec,lambda0,gamma,xvec_hole,yvec,zvec);

                gammad_vec(iterhole-na)=gammad;
                gammab_vec(iterhole-na)=gammab;
            end
            
            gammab_av(iternperp,iterL)=sum(2*abs(gvec(na+1:end)).^2.*gammab_vec);
            gammad_av(iternperp,iterL)=sum(2*abs(gvec(na+1:end)).^2.*gammad_vec);
            
        else 
            
            gammab_av(iternperp,iterL)=NaN;
            gammad_av(iternperp,iterL)=NaN;
            
        end
    end
end

delete(gcp)

save('hole_scan.mat','gammad_av','gammab_av')

%%

load('Compress.mat')
load('hole_scan.mat')

nperpvec=[4,8,12,16,20];
Nnperp=length(nperpvec);

xiopt=xiopt_cell{5};
omega0opt=omega0opt_cell{5};

Nxi=size(xiopt,2);

gvec_mean=zeros(Nnperp,Nxi);
gammabvec=zeros(Nnperp,Nxi);

for iternperp=1:Nnperp
    
    nperp=nperpvec(iternperp);
    na=nperp^2;
    
    for iterL=1:Nxi

        if 0.5*xiopt(iternperp,iterL)*(nperp*dperp)^2 > 1

            L=0.5*xiopt(iternperp,iterL)*(nperp*dperp)^2/lambda0+0.1;

            omega0=omega0opt(iternperp,iterL);

            zR=pi*omega0^2/lambda0;

            xvec=setXpos(nperp,dperp);
            yvec=setYpos(nperp,dperp);
            zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

            gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
            
            gvec_mean(iternperp,iterL)=sum(abs(gvec).^4);

            [gammad,gammabvec(iternperp,iterL)]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

        else 
            
            gvec_mean(iternperp,iterL)=NaN;
            gammabvec(iternperp,iterL)=NaN;
            
        end
    end
end

%% Plots


myBlue=[0.6,0.8,1];

figure
iterdperp=5;
ratiovec=ratiovec_cell{iterdperp};
xiopt=xiopt_cell{iterdperp};

for iternperp=1:Nnperp
    xi=xiopt(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    ratio=ratiovec(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    ratio_hole=gammad_av(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1)./gammab_av(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);

    h=plot(xi,ratio,'.-');
    h.MarkerSize=25;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    hold on

    h=plot(xi,ratio_hole,'.-');
    h.MarkerSize=25;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    hold on

    vmean=gvec_mean(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    gb=gammabvec(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    gd=ratio.*gb;
    gd_th=gd.*(1-2*vmean)+vmean*gamma;
    gb_th=gb.*(1-2*vmean)+vmean*gamma;

    h=plot(xi,gd_th./gb_th,'--');
    h.MarkerSize=25;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    hold on
end
chH = get(gca,'Children');
set(gca,'Children',[chH(end);chH(1:end-1)])

set(gca,'FontSize',22)
ax = gca;
ax.XLim = [1e-2, 1e1];
ax.XTick=[1e-2,1e-1,1e0,1e1];
ax.YTick=[1e-5,1e-4,1e-3,1e-2,1e-1,1e0];
ax.XScale='log';
ax.YScale='log';

ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
grid on
xlabel('$L\lambda_0/L_\perp^2$')
ylabel('$\langle\gamma_d^{(1 hole)}\rangle/\langle\gamma_b^{(1 hole)}\rangle$')
   
%%
myBlue=[0.6,0.8,1];

figure
iterdperp=5;
ratiovec=ratiovec_cell{iterdperp};
xiopt=xiopt_cell{iterdperp};

nperpvec=[4,8,12,16,20];
Nnperp=length(nperpvec);

for iternperp=1:Nnperp
    xi=xiopt(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    ratio=ratiovec(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    ratio_hole=gammad_av(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1)./gammab_av(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);

    h=plot(xi,nperpvec(iternperp)^2*(ratio_hole-ratio).*sqrt(xi),'.-');
    h.MarkerSize=25;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    hold on

end
chH = get(gca,'Children');
set(gca,'Children',[chH(end);chH(1:end-1)])

set(gca,'FontSize',22)
ax = gca;
ax.XLim = [1e-2, 1e1];
ax.XTick=[1e-2,1e-1,1e0,1e1];
ax.YLim=[1e-1,1e1];
ax.XScale='log';
ax.YScale='log';

ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
grid on
xlabel('$L\lambda_0/L_\perp^2$')
ylabel('$\sqrt{L\lambda_0/L_\perp^2}n_\perp^2 (\langle\gamma_d^{(1 hole)}\rangle/\langle\gamma_b^{(1 hole)}\rangle-\gamma_d/\gamma_b)$')
   