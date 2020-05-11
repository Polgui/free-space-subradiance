%% gammad/gammab
clearvars


lambda0=1;


nperp=15;
dperp=0.8;
L=15;
gamma=1e1;
 
na=nperp^2;

omega0=find_omega0opt(lambda0,nperp,dperp,L);
zR=pi*omega0^2/lambda0;


xvec=setXpos(nperp,dperp);
yvec=setYpos(nperp,dperp);
zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec,0*lambda0/4);

gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
[gammad,gammab,Deltalight,overlaps,Deltabright,ueig,ueig2,lam,u]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

gammad;
gammab;
ratiovec=gammad/gammab

%% Compute Efield

xmin=-10;
xmax=-xmin;
dx=0.05;

xgrid=(xmin:dx:xmax);

zmin=-max(15,1.5*L);
zmax=-zmin;
dz=0.01;

zgrid=(zmin:dz:zmax);

[X,Z]=meshgrid(xgrid,zgrid);
Y=zeros(size(X));

Efield_x=zeros(size(X));
Efield_y=zeros(size(X));
Efield_z=zeros(size(X));

Efield2_x=zeros(size(X));
Efield2_y=zeros(size(X));
Efield2_z=zeros(size(X));

Rat=0*dperp/4;

for itera=1:2*na
    Rmat=(X-xvec(itera)).^2+(Y-yvec(itera)).^2+(Z-zvec(itera)).^2;
    
    Efield_x=Efield_x+ueig(itera)*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(Rmat>Rat);
    Efield_y=Efield_y+ueig(itera)*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(Rmat>Rat);
    Efield_z=Efield_z+ueig(itera)*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(Rmat>Rat);
    
    Efield2_x=Efield2_x+ueig2(itera)*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(Rmat>Rat);
    Efield2_y=Efield2_y+ueig2(itera)*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(Rmat>Rat);
    Efield2_z=Efield2_z+ueig2(itera)*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(Rmat>Rat);
end

Efield=sqrt(Efield_x.^2+Efield_y.^2+Efield_z.^2);
Efield2=sqrt(Efield2_x.^2+Efield2_y.^2+Efield2_z.^2);

%%

figure('Position',[10,10,1500,10+1.2*(1200-10)*(xmax-xmin)/(zmax-zmin)])
set(gca,'color','none')
ax1=gca;

threshold=0.;
maxfield=max(max(abs(Efield)));
maxfield2=max(max(abs(Efield2)));
Emax=max(maxfield,maxfield2);
Enorm=abs(Efield/Emax);
Enorm=Enorm.*(Enorm>threshold);
h=imagesc(ax1,zgrid+L,xgrid,Enorm');
alpha('0.75');
hold on
colormap(ax1,(brewermap([],'Blues')))
caxis(ax1,[0,1]);
box off
ax1.XTick=[];
ax1.YTick=[];
axis equal

ax2=axes;
ax2.FontName = 'LaTeX';
ax2.Title.Interpreter = 'LaTeX';
ax2.XLabel.Interpreter = 'LaTeX';
ax2.YLabel.Interpreter = 'LaTeX';
ax2.TickLabelInterpreter = 'LaTeX';

ueig=ueig/max(abs(ueig));
% for itera=nperp*nperp/2+1:nperp*(nperp/2+1)
%     scatter(ax2,zvec(itera)+L,xvec(itera),300,abs(ueig(itera))^2,'filled','MarkerEdgeColor',0.9*[1,1,1],'LineWidth',2);
%     hold on
%     scatter(ax2,zvec(itera+na)+L,xvec(itera+na),300,abs(ueig(itera))^2,'filled','MarkerEdgeColor',0.9*[1,1,1],'LineWidth',2);
% end
axis equal

linkaxes([ax1,ax2])

colormap(ax2,brewermap([],'Reds'))
caxis(ax2,[min(abs(ueig).^2)-0.1*max(abs(ueig).^2),1.25*max(abs(ueig).^2)]);

box off
set([ax1,ax2],'FontSize',36);
set([ax1,ax2],'FontName','LaTeX');
set([ax1,ax2],'XTick',[]);
set([ax1,ax2],'YTick',[]);
set([ax1,ax2],'Box','off');
ax2.Visible='off';
ax1.TickLabelInterpreter = 'LaTeX';
axis equal

%%

figure('Position',[10,10,1500,10+1.2*(1200-10)*(xmax-xmin)/(zmax-zmin)])
set(gca,'color','none')
ax1=gca;

Enorm=abs(Efield2/Emax);
Enorm=Enorm.*(Enorm>threshold);
h=imagesc(ax1,zgrid+L,xgrid,Enorm');
hold on
colormap(ax1,(brewermap([],'Blues')))
box off
ax1.XTick=[];
ax1.YTick=[];
axis equal
alpha('0.75');

ax2=axes;
ax2.FontName = 'LaTeX';
ax2.Title.Interpreter = 'LaTeX';
ax2.XLabel.Interpreter = 'LaTeX';
ax2.YLabel.Interpreter = 'LaTeX';
ax2.TickLabelInterpreter = 'LaTeX';

ueig2=ueig2/max(abs(ueig2));
thresholdedge=0.7;
for itera=nperp*nperp/2+1:nperp*(nperp/2+1)
    scatter(ax2,zvec(itera)+L,xvec(itera),30,abs(ueig2(itera))^2,'filled','MarkerEdgeColor',0.9*[1,1,1],'LineWidth',2);
    hold on
    scatter(ax2,zvec(itera+na)+L,xvec(itera+na),300,abs(ueig2(itera))^2,'filled','MarkerEdgeColor',0.9*[1,1,1],'LineWidth',2);
end
axis equal


linkaxes([ax1,ax2])

colormap(ax2,brewermap([],'Reds'))
caxis(ax2,[min(abs(ueig2).^2)-0.1*max(abs(ueig).^2),1.25*max(abs(ueig2).^2)]);

box off
set([ax1,ax2],'FontSize',36);
set([ax1,ax2],'FontName','LaTeX');
set([ax1,ax2],'XTick',[]);
set([ax1,ax2],'YTick',[]);
set([ax1,ax2],'Box','off');
ax2.Visible='off';
ax1.TickLabelInterpreter = 'LaTeX';
axis equal

plot(ax2,zgrid+L,1.5*omega0*sqrt(1+(zgrid/zR).^2),'--r')
plot(ax2,zgrid+L,-1.5*omega0*sqrt(1+(zgrid/zR).^2),'--r')

%% Compute Efield hole

xmin=-10;
xmax=-xmin;
dx=0.1;

xgrid=(xmin:dx:xmax);

zmin=-1.5*L;
zmax=-zmin;
dz=0.01;

zgrid=(zmin:dz:zmax);

[X,Z]=meshgrid(xgrid,zgrid);
Y=zeros(size(X));

Efield_x=zeros(size(X));
Efield_y=zeros(size(X));
Efield_z=zeros(size(X));

Efield2_x=zeros(size(X));
Efield2_y=zeros(size(X));
Efield2_z=zeros(size(X));

Rat=dperp/4;

for itera=1:2*na
    Rmat=(X-xvec(itera)).^2+(Y-yvec(itera)).^2+(Z-zvec(itera)).^2;
    
    Efield_x=Efield_x+ueig(itera)*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(itera~=na+nperp*nperp/2+nperp/2);
    Efield_y=Efield_y+ueig(itera)*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(itera~=na+nperp*nperp/2+nperp/2);
    Efield_z=Efield_z+ueig(itera)*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(itera~=na+nperp*nperp/2+nperp/2);
    
    Efield2_x=Efield2_x+ueig2(itera)*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(itera~=na+nperp*nperp/2+nperp/2);
    Efield2_y=Efield2_y+ueig2(itera)*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(itera~=na+nperp*nperp/2+nperp/2);
    Efield2_z=Efield2_z+ueig2(itera)*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(itera~=na+nperp*nperp/2+nperp/2);
end

Efield=sqrt(Efield_x.^2+Efield_y.^2+Efield_z.^2);
Efield2=sqrt(Efield2_x.^2+Efield2_y.^2+Efield2_z.^2);
%%
figure('Position',[10,10,1500,10+1.2*(1200-10)*(xmax-xmin)/(zmax-zmin)])
set(gca,'color','none')
ax1=gca;

threshold=0.;
maxfield=max(max(abs(Efield)));
Enorm=abs(Efield/maxfield);
Enorm=Enorm.*(Enorm>threshold);
h=imagesc(ax1,zgrid+L,xgrid,(Enorm').^(1));
hold on
colormap(ax1,(brewermap([],'Blues')))
box off
ax1.XTick=[];
ax1.YTick=[];
axis equal

ax2=axes;
ax2.FontName = 'LaTeX';
ax2.Title.Interpreter = 'LaTeX';
ax2.XLabel.Interpreter = 'LaTeX';
ax2.YLabel.Interpreter = 'LaTeX';
ax2.TickLabelInterpreter = 'LaTeX';

for itera=nperp*nperp/2+1:nperp*(nperp/2+1)
    scatter(ax2,zvec(itera)+L,-xvec(itera),120,abs(ueig(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
    hold on
    if itera~=nperp*nperp/2+nperp/2
        scatter(ax2,zvec(itera+na)+L,-xvec(itera+na),120,abs(ueig(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
    else
        %scatter(ax2,zvec(itera+na)+L,-xvec(itera+na),120,'filled','MarkerEdgeColor','k','MarkerFaceColor','w');
    end
end
axis equal

linkaxes([ax1,ax2])

colormap(ax2,brewermap([],'Reds'))
caxis(ax2,[min(abs(ueig).^2)-0.1*max(abs(ueig).^2),1.25*max(abs(ueig).^2)]);

box off
set([ax1,ax2],'FontSize',36);
set([ax1,ax2],'FontName','LaTeX');
set([ax1,ax2],'XTick',[]);
set([ax1,ax2],'YTick',[]);
set([ax1,ax2],'Box','off');
ax2.Visible='off';
ax1.TickLabelInterpreter = 'LaTeX';
axis equal

%%
figure('Position',[10,10,1500,10+1.2*(1200-10)*(xmax-xmin)/(zmax-zmin)])
set(gca,'color','none')
ax1=gca;

maxfield=max(max(abs(Efield2)));
Enorm=abs(Efield2/maxfield);
h=imagesc(ax1,zgrid+L,xgrid,Enorm');
hold on
colormap(ax1,(brewermap([],'Blues')))
box off
ax1.XTick=[];
ax1.YTick=[];
axis equal

ax2=axes;
ax2.FontName = 'LaTeX';
ax2.Title.Interpreter = 'LaTeX';
ax2.XLabel.Interpreter = 'LaTeX';
ax2.YLabel.Interpreter = 'LaTeX';
ax2.TickLabelInterpreter = 'LaTeX';

for itera=nperp*nperp/2+1:nperp*(nperp/2+1)
    scatter(ax2,zvec(itera)+L,xvec(itera),120,abs(ueig2(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
    hold on
    if itera~=nperp*nperp/2+nperp/2
        scatter(ax2,zvec(itera+na)+L,-xvec(itera+na),120,abs(ueig2(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
    else
        %scatter(ax2,zvec(itera+na)+L,-xvec(itera+na),120,'filled','MarkerEdgeColor','k','MarkerFaceColor','w');
    end
end
axis equal

linkaxes([ax1,ax2])

colormap(ax2,brewermap([],'Reds'))
caxis(ax2,[min(abs(ueig2).^2)-0.1*max(abs(ueig).^2),1.25*max(abs(ueig2).^2)]);

box off
set([ax1,ax2],'FontSize',36);
set([ax1,ax2],'FontName','LaTeX');
set([ax1,ax2],'XTick',[]);
set([ax1,ax2],'YTick',[]);
set([ax1,ax2],'Box','off');
ax2.Visible='off';
ax1.TickLabelInterpreter = 'LaTeX';
axis equal


%% Compute Efield 1 layer

xmin=-10;
xmax=-xmin;
dx=0.1;

xgrid=(xmin:dx:xmax);

zmin=-L-15;
zmax=-zmin;
dz=0.01;

zgrid=(zmin:dz:zmax);

[X,Z]=meshgrid(xgrid,zgrid);
Y=zeros(size(X));

Efield_x=zeros(size(X));
Efield_y=zeros(size(X));
Efield_z=zeros(size(X));


Rat=0*dperp/4;

for itera=1:na
    Rmat=(X-xvec(itera)).^2+(Y-yvec(itera)).^2+(Z-zvec(itera)).^2;
    
    Efield_x=Efield_x+ueig(itera)*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(Rmat>Rat);
    Efield_y=Efield_y+ueig(itera)*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(Rmat>Rat);
    Efield_z=Efield_z+ueig(itera)*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(Rmat>Rat);
    
end

Efield=sqrt(Efield_x.^2+Efield_y.^2+Efield_z.^2);

figure('Position',[10,10,1500,10+1.2*(1200-10)*(xmax-xmin)/(zmax-zmin)])
set(gca,'color','none')
ax1=gca;

maxfield=max(max(abs(Efield)));
Enorm=abs(Efield/maxfield);
h=imagesc(ax1,zgrid+L,xgrid,Enorm');
hold on
colormap(ax1,(brewermap([],'Blues')))
box off
ax1.XTick=[];
ax1.YTick=[];
axis equal

ax2=axes;
ax2.FontName = 'LaTeX';
ax2.Title.Interpreter = 'LaTeX';
ax2.XLabel.Interpreter = 'LaTeX';
ax2.YLabel.Interpreter = 'LaTeX';
ax2.TickLabelInterpreter = 'LaTeX';

for itera=nperp*nperp/2+1:nperp*(nperp/2+1)
    scatter(ax2,zvec(itera)+L,xvec(itera),120,abs(ueig(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
    hold on
end
axis equal

linkaxes([ax1,ax2])

colormap(ax2,brewermap([],'Reds'))
caxis(ax2,[min(abs(ueig).^2)-0.1*max(abs(ueig).^2),1.25*max(abs(ueig).^2)]);

box off
set([ax1,ax2],'FontSize',36);
set([ax1,ax2],'FontName','LaTeX');
set([ax1,ax2],'XTick',[]);
set([ax1,ax2],'YTick',[]);
set([ax1,ax2],'Box','off');
ax2.Visible='off';
ax1.TickLabelInterpreter = 'LaTeX';
axis equal

%% Compute Efield 1 layer flat no phase

xmin=-10;
xmax=-xmin;
dx=0.1;

xgrid=(xmin:dx:xmax);

zmin=-L-15;
zmax=-zmin;
dz=0.01;

zgrid=(zmin:dz:zmax);

[X,Z]=meshgrid(xgrid,zgrid);
Y=zeros(size(X));

Efield_x=zeros(size(X));
Efield_y=zeros(size(X));
Efield_z=zeros(size(X));


Rat=dperp/4;
zat=-L;

for itera=1:na
    Rmat=(X-xvec(itera)).^2+(Y-yvec(itera)).^2+(Z-zat).^2;
    
    Efield_x=Efield_x+ueig(itera)*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
    Efield_y=Efield_y+ueig(itera)*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
    Efield_z=Efield_z+ueig(itera)*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
    
end

Efield=sqrt(Efield_x.^2+Efield_y.^2+Efield_z.^2);

%%
figure('Position',[10,10,1500,10+1.2*(1200-10)*(xmax-xmin)/(zmax-zmin)])
set(gca,'color','none')
ax1=gca;

maxfield=max(max(abs(Efield)));
Enorm=abs(Efield/maxfield);
h=imagesc(ax1,zgrid+L,xgrid,Enorm');
hold on
colormap(ax1,(brewermap([],'Blues')))
box off
ax1.XTick=[];
ax1.YTick=[];
axis equal

ax2=axes;
ax2.FontName = 'LaTeX';
ax2.Title.Interpreter = 'LaTeX';
ax2.XLabel.Interpreter = 'LaTeX';
ax2.YLabel.Interpreter = 'LaTeX';
ax2.TickLabelInterpreter = 'LaTeX';

for itera=nperp*nperp/2+1:nperp*(nperp/2+1)
    scatter(ax2,zat+L,xvec(itera),120,abs(ueig(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
    hold on
end
axis equal

linkaxes([ax1,ax2])

colormap(ax2,brewermap([],'Reds'))
caxis(ax2,[min(abs(ueig).^2)-0.1*max(abs(ueig).^2),1.25*max(abs(ueig).^2)]);

box off
set([ax1,ax2],'FontSize',36);
set([ax1,ax2],'FontName','LaTeX');
set([ax1,ax2],'XTick',[]);
set([ax1,ax2],'YTick',[]);
set([ax1,ax2],'Box','off');
ax2.Visible='off';
ax1.TickLabelInterpreter = 'LaTeX';
axis equal

%% Compute Efield 1 layer flat gaussian

xmin=-10;
xmax=-xmin;
dx=0.1;

xgrid=(xmin:dx:xmax);

zmin=-L-15;
zmax=-zmin;
dz=0.01;

zgrid=(zmin:dz:zmax);

[X,Z]=meshgrid(xgrid,zgrid);
Y=zeros(size(X));

Efield_x=zeros(size(X));
Efield_y=zeros(size(X));
Efield_z=zeros(size(X));

Rat=dperp/4;
zat=L;

for itera=na+1:2*na
    Rmat=(X-xvec(itera)).^2+(Y-yvec(itera)).^2+(Z-zat).^2;
    
    r2=xvec(itera)^2+yvec(itera)^2;
    zR=pi*omega0^2/lambda0;
    omegaz=omega0*sqrt(1+(zat/zR).^2);
    Rz=zat*(1+(zat./L).^2);
    psiz=atan(zat/zR);

    GaussianBeam=sqrt(2/pi)*(1./omegaz).*exp(-r2./(omegaz.^2)).*exp(-1j*psiz).*exp(2*pi*1j*(zat+r2./(2*Rz))/lambda0);


    
    Efield_x=Efield_x+GaussianBeam*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
    Efield_y=Efield_y+GaussianBeam*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
    Efield_z=Efield_z+GaussianBeam*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
   
end

Efield=sqrt(Efield_x.^2+Efield_y.^2+Efield_z.^2);

figure('Position',[10,10,1500,10+1.2*(1200-10)*(xmax-xmin)/(zmax-zmin)])
set(gca,'color','none')
ax1=gca;

maxfield=max(max(abs(Efield)));
Enorm=abs(Efield/maxfield);
h=imagesc(ax1,zgrid+L,xgrid,Enorm');
hold on
colormap(ax1,(brewermap([],'Blues')))
box off
ax1.XTick=[];
ax1.YTick=[];
axis equal

ax2=axes;
ax2.FontName = 'LaTeX';
ax2.Title.Interpreter = 'LaTeX';
ax2.XLabel.Interpreter = 'LaTeX';
ax2.YLabel.Interpreter = 'LaTeX';
ax2.TickLabelInterpreter = 'LaTeX';

for itera=nperp*nperp/2+1:nperp*(nperp/2+1)
    scatter(ax2,zat+L,xvec(itera),120,abs(ueig(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
    hold on
end
axis equal

linkaxes([ax1,ax2])

colormap(ax2,brewermap([],'Reds'))
caxis(ax2,[min(abs(ueig).^2)-0.1*max(abs(ueig).^2),1.25*max(abs(ueig).^2)]);

box off
set([ax1,ax2],'FontSize',36);
set([ax1,ax2],'FontName','LaTeX');
set([ax1,ax2],'XTick',[]);
set([ax1,ax2],'YTick',[]);
set([ax1,ax2],'Box','off');
ax2.Visible='off';
ax1.TickLabelInterpreter = 'LaTeX';
axis equal


%% Compute Efield 2 layers flat no phase

xmin=-10;
xmax=-xmin;
dx=0.1;

xgrid=(xmin:dx:xmax);

zmin=-L-15;
zmax=-zmin;
dz=0.01;

zgrid=(zmin:dz:zmax);

[X,Z]=meshgrid(xgrid,zgrid);
Y=zeros(size(X));

Efield_x=zeros(size(X));
Efield_y=zeros(size(X));
Efield_z=zeros(size(X));

Efield2_x=zeros(size(X));
Efield2_y=zeros(size(X));
Efield2_z=zeros(size(X));

Rat=dperp/4;
zat=-L;

for itera=1:na
    Rmat=(X-xvec(itera)).^2+(Y-yvec(itera)).^2+(Z-zat).^2;
    
    Efield_x=Efield_x+ueig(itera)*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
    Efield_y=Efield_y+ueig(itera)*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
    Efield_z=Efield_z+ueig(itera)*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
    
    Efield2_x=Efield2_x+ueig2(itera)*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
    Efield2_y=Efield2_y+ueig2(itera)*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
    Efield2_z=Efield2_z+ueig2(itera)*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
end
for itera=na+1:2*na
    Rmat=(X-xvec(itera)).^2+(Y-yvec(itera)).^2+(Z+zat).^2;
    
    Efield_x=Efield_x+ueig(itera)*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z+zat).*(Rmat>Rat);
    Efield_y=Efield_y+ueig(itera)*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z+zat).*(Rmat>Rat);
    Efield_z=Efield_z+ueig(itera)*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z+zat).*(Rmat>Rat);
    
    Efield2_x=Efield2_x+ueig2(itera)*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z+zat).*(Rmat>Rat);
    Efield2_y=Efield2_y+ueig2(itera)*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z+zat).*(Rmat>Rat);
    Efield2_z=Efield2_z+ueig2(itera)*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z+zat).*(Rmat>Rat);
end

Efield=sqrt(Efield_x.^2+Efield_y.^2+Efield_z.^2);
Efield2=sqrt(Efield2_x.^2+Efield2_y.^2+Efield2_z.^2);

figure('Position',[10,10,1500,10+1.2*(1200-10)*(xmax-xmin)/(zmax-zmin)])
set(gca,'color','none')
ax1=gca;

maxfield=max(max(abs(Efield2)));
Enorm=abs(Efield2/maxfield);
h=imagesc(ax1,zgrid+L,xgrid,Enorm');
hold on
colormap(ax1,(brewermap([],'Blues')))
box off
ax1.XTick=[];
ax1.YTick=[];
axis equal

ax2=axes;
ax2.FontName = 'LaTeX';
ax2.Title.Interpreter = 'LaTeX';
ax2.XLabel.Interpreter = 'LaTeX';
ax2.YLabel.Interpreter = 'LaTeX';
ax2.TickLabelInterpreter = 'LaTeX';

for itera=nperp*nperp/2+1:nperp*(nperp/2+1)
    scatter(ax2,zat+L,xvec(itera),120,abs(ueig(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
    hold on
    scatter(ax2,-zat+L,xvec(na+itera),120,abs(ueig(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
end
axis equal

linkaxes([ax1,ax2])

colormap(ax2,brewermap([],'Reds'))
caxis(ax2,[min(abs(ueig).^2)-0.1*max(abs(ueig).^2),1.25*max(abs(ueig).^2)]);

box off
set([ax1,ax2],'FontSize',36);
set([ax1,ax2],'FontName','LaTeX');
set([ax1,ax2],'XTick',[]);
set([ax1,ax2],'YTick',[]);
set([ax1,ax2],'Box','off');
ax2.Visible='off';
ax1.TickLabelInterpreter = 'LaTeX';
axis equal

%% Compute Efield 2 layers flat gaussian

xmin=-10;
xmax=-xmin;
dx=0.1;

xgrid=(xmin:dx:xmax);

zmin=-L-15;
zmax=-zmin;
dz=0.01;

zgrid=(zmin:dz:zmax);

[X,Z]=meshgrid(xgrid,zgrid);
Y=zeros(size(X));

Efield_x=zeros(size(X));
Efield_y=zeros(size(X));
Efield_z=zeros(size(X));


Rat=dperp/4;
zat=-L;

for itera=1:na
    Rmat=(X-xvec(itera)).^2+(Y-yvec(itera)).^2+(Z-zat).^2;
    
    r2=xvec(itera)^2+yvec(itera)^2;
    zR=pi*omega0^2/lambda0;
    omegaz=omega0*sqrt(1+(-L/zR).^2);
    Rz=-L.*(1+(-zR./L).^2);
    psiz=atan(-L/zR);

    GaussianBeam=sqrt(2/pi)*(1./omegaz).*exp(-r2./(omegaz.^2)).*exp(-1j*psiz).*exp(2*pi*1j*(-L+r2./(2*Rz))/lambda0);


    
    Efield_x=Efield_x+GaussianBeam*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
    Efield_y=Efield_y+GaussianBeam*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
    Efield_z=Efield_z+GaussianBeam*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z-zat).*(Rmat>Rat);
   
end

for itera=na+1:2*na
    Rmat=(X-xvec(itera)).^2+(Y-yvec(itera)).^2+(Z+zat).^2;
    
    r2=xvec(itera)^2+yvec(itera)^2;
    zR=pi*omega0^2/lambda0;
    omegaz=omega0*sqrt(1+(L/zR).^2);
    Rz=L.*(1+(zR./L).^2);
    psiz=atan(L/zR);

    GaussianBeam=sqrt(2/pi)*(1./omegaz).*exp(-r2./(omegaz.^2)).*exp(-1j*psiz).*exp(2*pi*1j*(L+r2./(2*Rz))/lambda0);


    
    Efield_x=Efield_x-GaussianBeam*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z+zat).*(Rmat>Rat);
    Efield_y=Efield_y-GaussianBeam*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z+zat).*(Rmat>Rat);
    Efield_z=Efield_z-GaussianBeam*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z+zat).*(Rmat>Rat);
   
end

Efield=sqrt(Efield_x.^2+Efield_y.^2+Efield_z.^2);

figure('Position',[10,10,1500,10+1.2*(1200-10)*(xmax-xmin)/(zmax-zmin)])
set(gca,'color','none')
ax1=gca;

maxfield=max(max(abs(Efield)));
Enorm=abs(Efield/maxfield);
h=imagesc(ax1,zgrid+L,xgrid,Enorm');
hold on
colormap(ax1,(brewermap([],'Blues')))
box off
ax1.XTick=[];
ax1.YTick=[];
axis equal

ax2=axes;
ax2.FontName = 'LaTeX';
ax2.Title.Interpreter = 'LaTeX';
ax2.XLabel.Interpreter = 'LaTeX';
ax2.YLabel.Interpreter = 'LaTeX';
ax2.TickLabelInterpreter = 'LaTeX';

for itera=nperp*nperp/2+1:nperp*(nperp/2+1)
    scatter(ax2,zat+L,xvec(itera),120,abs(ueig(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
    hold on
    scatter(ax2,-zat+L,xvec(na+itera),120,abs(ueig(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
end
axis equal

linkaxes([ax1,ax2])

colormap(ax2,brewermap([],'Reds'))
caxis(ax2,[min(abs(ueig).^2)-0.1*max(abs(ueig).^2),1.25*max(abs(ueig).^2)]);

box off
set([ax1,ax2],'FontSize',36);
set([ax1,ax2],'FontName','LaTeX');
set([ax1,ax2],'XTick',[]);
set([ax1,ax2],'YTick',[]);
set([ax1,ax2],'Box','off');
ax2.Visible='off';
ax1.TickLabelInterpreter = 'LaTeX';
axis equal

%% Compute Efield 2 holes

xmin=-10;
xmax=-xmin;
dx=0.1;

xgrid=(xmin:dx:xmax);

zmin=-1.5*L;
zmax=-zmin;
dz=0.01;

zgrid=(zmin:dz:zmax);

[X,Z]=meshgrid(xgrid,zgrid);
Y=zeros(size(X));

Efield_x=zeros(size(X));
Efield_y=zeros(size(X));
Efield_z=zeros(size(X));

Efield2_x=zeros(size(X));
Efield2_y=zeros(size(X));
Efield2_z=zeros(size(X));

Rat=dperp/4;

for itera=1:2*na
    Rmat=(X-xvec(itera)).^2+(Y-yvec(itera)).^2+(Z-zvec(itera)).^2;
    
    Efield_x=Efield_x+ueig(itera)*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(itera~=na+nperp*nperp/2+nperp/2).*(itera~=na+nperp*nperp/2+nperp/2+1);
    Efield_y=Efield_y+ueig(itera)*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(itera~=na+nperp*nperp/2+nperp/2).*(itera~=na+nperp*nperp/2+nperp/2+1);
    Efield_z=Efield_z+ueig(itera)*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(itera~=na+nperp*nperp/2+nperp/2).*(itera~=na+nperp*nperp/2+nperp/2+1);
    
    Efield2_x=Efield2_x+ueig2(itera)*GreensTensor_x(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(itera~=na+nperp*nperp/2+nperp/2).*(itera~=na+nperp*nperp/2+nperp/2+1);
    Efield2_y=Efield2_y+ueig2(itera)*GreensTensor_y(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(itera~=na+nperp*nperp/2+nperp/2).*(itera~=na+nperp*nperp/2+nperp/2+1);
    Efield2_z=Efield2_z+ueig2(itera)*GreensTensor_z(lambda0,X-xvec(itera),Y-yvec(itera),Z-zvec(itera)).*(itera~=na+nperp*nperp/2+nperp/2).*(itera~=na+nperp*nperp/2+nperp/2+1);
end

Efield=sqrt(Efield_x.^2+Efield_y.^2+Efield_z.^2);
Efield2=sqrt(Efield2_x.^2+Efield2_y.^2+Efield2_z.^2);
%%
figure('Position',[10,10,1500,10+1.2*(1200-10)*(xmax-xmin)/(zmax-zmin)])
set(gca,'color','none')
ax1=gca;

threshold=0.;
maxfield=max(max(abs(Efield)));
Enorm=abs(Efield/maxfield);
Enorm=Enorm.*(Enorm>threshold);
h=imagesc(ax1,zgrid+L,xgrid,Enorm');
hold on
colormap(ax1,(brewermap([],'Blues')))
box off
ax1.XTick=[];
ax1.YTick=[];
axis equal

ax2=axes;
ax2.FontName = 'LaTeX';
ax2.Title.Interpreter = 'LaTeX';
ax2.XLabel.Interpreter = 'LaTeX';
ax2.YLabel.Interpreter = 'LaTeX';
ax2.TickLabelInterpreter = 'LaTeX';

for itera=nperp*nperp/2+1:nperp*(nperp/2+1)
    scatter(ax2,zvec(itera)+L,-xvec(itera),120,abs(ueig(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
    hold on
    if itera~=nperp*nperp/2+nperp/2 && itera~=nperp*nperp/2+nperp/2+1
        scatter(ax2,zvec(itera+na)+L,-xvec(itera+na),120,abs(ueig(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
    else
        %scatter(ax2,zvec(itera+na)+L,-xvec(itera+na),120,'filled','MarkerEdgeColor','k','MarkerFaceColor','w');
    end
end
axis equal

linkaxes([ax1,ax2])

colormap(ax2,brewermap([],'Reds'))
caxis(ax2,[min(abs(ueig).^2)-0.1*max(abs(ueig).^2),1.25*max(abs(ueig).^2)]);

box off
set([ax1,ax2],'FontSize',36);
set([ax1,ax2],'FontName','LaTeX');
set([ax1,ax2],'XTick',[]);
set([ax1,ax2],'YTick',[]);
set([ax1,ax2],'Box','off');
ax2.Visible='off';
ax1.TickLabelInterpreter = 'LaTeX';
axis equal

%%
figure('Position',[10,10,1500,10+1.2*(1200-10)*(xmax-xmin)/(zmax-zmin)])
set(gca,'color','none')
ax1=gca;

maxfield=max(max(abs(Efield2)));
Enorm=abs(Efield2/maxfield);
h=imagesc(ax1,zgrid+L,xgrid,Enorm');
hold on
colormap(ax1,(brewermap([],'Blues')))
box off
ax1.XTick=[];
ax1.YTick=[];
axis equal

ax2=axes;
ax2.FontName = 'LaTeX';
ax2.Title.Interpreter = 'LaTeX';
ax2.XLabel.Interpreter = 'LaTeX';
ax2.YLabel.Interpreter = 'LaTeX';
ax2.TickLabelInterpreter = 'LaTeX';

for itera=nperp*nperp/2+1:nperp*(nperp/2+1)
    scatter(ax2,zvec(itera)+L,xvec(itera),120,abs(ueig2(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
    hold on
    if itera~=nperp*nperp/2+nperp/2
        scatter(ax2,zvec(itera+na)+L,-xvec(itera+na),120,abs(ueig2(itera))^2,'filled','MarkerEdgeColor',[1,1,1]);
    else
        %scatter(ax2,zvec(itera+na)+L,-xvec(itera+na),120,'filled','MarkerEdgeColor','k','MarkerFaceColor','w');
    end
end
axis equal

linkaxes([ax1,ax2])

colormap(ax2,brewermap([],'Reds'))
caxis(ax2,[min(abs(ueig2).^2)-0.1*max(abs(ueig).^2),1.25*max(abs(ueig2).^2)]);

box off
set([ax1,ax2],'FontSize',36);
set([ax1,ax2],'FontName','LaTeX');
set([ax1,ax2],'XTick',[]);
set([ax1,ax2],'YTick',[]);
set([ax1,ax2],'Box','off');
ax2.Visible='off';
ax1.TickLabelInterpreter = 'LaTeX';
axis equal
