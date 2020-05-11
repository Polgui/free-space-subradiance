clearvars
addpath('expmv-master')

nperpvec=[4,5,6,7,8];
Nnperp=length(nperpvec);

ximin=-1;
ximax=1;
Nxi=24;
xivec=10.^linspace(ximin,ximax,Nxi);

mycluster=parcluster();
mypool=parpool(mycluster);  

xioptmat=zeros(Nnperp,Nxi);
ratiomat=zeros(Nnperp,Nxi);
ratiomat2=zeros(Nnperp,Nxi);
Deltamat=zeros(Nnperp,Nxi);
Deltamat2=zeros(Nnperp,Nxi);
overlapmin=zeros(Nnperp,Nxi);
overlapmax=zeros(Nnperp,Nxi);

for iternperp=1:Nnperp
    disp(nperpvec(iternperp))
    parfor iterxi=1:Nxi
        nperp=nperpvec(iternperp);
        dperp=0.5;
        L=0.5*xivec(iterxi)*(nperp*dperp)^2;
        if L>1
            lambda0=1;
            gamma=1e1;
            omega0=find_omega0opt(lambda0,nperp,dperp,L);

            na=nperp^2;
            zR=pi*omega0^2/lambda0;

            xvec=setXpos(nperp,dperp);
            yvec=setYpos(nperp,dperp);
            zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

            gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);

            [gammad,gammab,Deltalight]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);
            [gammad2,eigvec2,sortoverlaps,Deltalight2]=getShift2(gvec,lambda0,gamma,xvec,yvec,zvec);
            
            ratiomat(iternperp,iterxi)=gammad/gammab;
            ratiomat2(iternperp,iterxi)=gammad2/(2*gammab);
            Deltamat(iternperp,iterxi)=Deltalight/gammab;
            Deltamat2(iternperp,iterxi)=Deltalight2/(2*gammab);
            xioptmat(iternperp,iterxi)=2*max(zvec)*lambda0/(nperp*dperp)^2;
            overlapmin(iternperp,iterxi)=sum(abs(sortoverlaps(1:end-1)).^2);
            overlapmax(iternperp,iterxi)=abs(sortoverlaps(end)).^2;
        else
            ratiomat(iternperp,iterxi)=NaN;
            ratiomat2(iternperp,iterxi)=NaN;
            Deltamat(iternperp,iterxi)=NaN;
            Deltamat2(iternperp,iterxi)=NaN;
            xioptmat(iternperp,iterxi)=2*L/(nperp*dperp)^2;
            overlapmin(iternperp,iterxi)=NaN;
            overlapmax(iternperp,iterxi)=NaN;
        end
    end
end

delete(gcp)

save('2exc_scan.mat','ratiomat','ratiomat2','Deltamat','Deltamat2','xioptmat','overlapmax','overlapmin')

%%

nperpvec=[4,5,6,7,8];
Nnperp=length(nperpvec);

ximin=-1;
ximax=1;
Nxi=24;
xivec=10.^linspace(ximin,ximax,Nxi);

gvec_mean=zeros(Nnperp,Nxi);
gammabvec=zeros(Nnperp,Nxi);
xiopt_mean=zeros(Nnperp,Nxi);

for iternperp=1:Nnperp
    
    dperp=0.5;
    gamma=10;
    lambda0=1;
    nperp=nperpvec(iternperp)
    na=nperp^2;
    
    for iterL=1:Nxi

        L=0.5*xivec(iterL)*(nperp*dperp)^2/lambda0+0.1;

        omega0=find_omega0opt(lambda0,nperp,dperp,L);

        zR=pi*omega0^2/lambda0;

        xvec=setXpos(nperp,dperp);
        yvec=setYpos(nperp,dperp);
        zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

        gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);

        gvec_mean(iternperp,iterL)=sum(abs(gvec).^4);

        [gammad,gammabvec(iternperp,iterL)]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);
        xiopt_mean(iternperp,iterL)=2*max(zvec)/(dperp*nperp)^2;
    end
end

save('gvec_mean.mat','gvec_mean','gammabvec','xiopt_mean')

%%

load('2exc_scan.mat')
load('gvec_mean.mat')

nperpvec=[4,5,6,7,8];
Nnperp=length(nperpvec);

figure

myBlue=[0.6,1,0.8];

mylegend=cell(Nnperp,1);
plot([0,1],[0,1],'w');
hold on
% for iternperp=1:Nnperp
%     xi=xioptmat(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
%     ratio=ratiomat(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
%     h=plot(xi,ratio,'.-');
%     h.MarkerSize=25;
%     h.LineWidth=2;
%     h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
%     hold on
% end

set(gca,'FontSize',22)
ax = gca;
ax.XLim = [1e-1, 1e1];
ax.XTick=[1e-1,1e0,1e1];
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
ylabel('$\gamma_d/\gamma_b$')


for iternperp=1:Nnperp
    xi=xioptmat(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    ratio=ratiomat2(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    h=plot(xi,ratio,'.-');
    h.MarkerSize=25;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-iternperp+2)/(Nnperp+1);
    hold on
end

for iternperp=1:Nnperp
    xi=xiopt_mean(iternperp,0.5*xiopt_mean(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    
    ratio=ratiomat(iternperp,0.5*xiopt_mean(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    vmean=gvec_mean(iternperp,0.5*xiopt_mean(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    gb=gammabvec(iternperp,0.5*xiopt_mean(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    gd=ratio.*gb;
    gd_th=gd.*(1-2*vmean)+vmean*gamma;
    gb_th=gb.*(1-2*vmean)+vmean*gamma;

    h=plot(xi,(gd_th)./(gb_th),'--');
    h.MarkerSize=25;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
end

%%

myBlue=[0.6,0.8,1];

Nxi=24;
figure
for iternperp=1:Nnperp
    for iterxi=1:Nxi
        h=plot(xioptmat(iternperp,iterxi),overlapmin(iternperp,iterxi),'.-');    
        h.MarkerSize=25;
        h.LineWidth=2;
        h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
        hold on
        h=plot(xioptmat(iternperp,iterxi),overlapmax(iternperp,iterxi),'.-');   
        h.MarkerSize=25;
        h.LineWidth=2;
        h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    end
end
set(gca,'FontSize',22)
ax = gca;
ax.XLim = [1e-1, 1e1];
ax.XTick=[1e-1,1e0,1e1];
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
ylabel('$O$')
