%% gammad/gammab
clearvars


lambda0=1;


nperp=20;
dperp=0.5;
L=5;
gamma=1e1;
 
na=nperp^2;

omega0=find_omega0opt(lambda0,nperp,dperp,L);
zR=pi*omega0^2/lambda0;


xvec=setXpos(nperp,dperp);
yvec=setYpos(nperp,dperp);
zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec,lambda0/4);

gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
[gammad,gammab,Deltalight,overlaps,Deltabright,ueig,ueig2,lam,u]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

gammad;
gammab;
ratiovec=gammad/gammab

%% Compute Efield

xmin=-8;
xmax=-xmin;
dx=0.04;

xgrid=(xmin:dx:xmax);

zmin=-1.5*L;
zmax=-zmin;
dz=0.01;

zgrid=(zmin:dz:zmax);

[X,Z]=meshgrid(xgrid,zgrid);
Y=zeros(size(X));

Efield=zeros(size(X));

for itera=1:na
    Efield=Efield+GreensTensor(lambda0,sqrt((X-xvec(itera)).^2+(Y-yvec(itera)).^2),Z-zvec(itera));
end
for itera=na+1:2*na
    Efield=Efield+GreensTensor(lambda0,sqrt((X-xvec(itera)).^2+(Y-yvec(itera)).^2),Z-zvec(itera));
end


%% Plot Efield
figure
set(gca,'color','none')

a=0.18;
l = [-.25 -.433 1.5];
for j=1:length(zvec)
  % Generate a sphere consisting of 20by 20 faces
  [x,y,z]=sphere;
  c = max(0,x*l(1) + y*l(2) + z*l(3));
  % use surf function to plot
  hSurface=surf(a*x+zvec(j)+L,a*y+xvec(j),a*z+yvec(j),c);
  hold on
  set(hSurface,'FaceColor','r', ...
      'FaceAlpha',0.5,'EdgeColor','none')
  daspect([1 1 1]);
end
%camlight
%lighting phong


hAxes = gca;
threshEfield=0.;
maxfield=max(max(abs(Efield)));
Enorm=abs(Efield/maxfield);
Efilter=Enorm.*(Enorm>threshEfield);
Efinal=(Enorm-threshEfield)/max(max(Enorm-threshEfield));
h=pcolor(Z+L,X,Efilter.^2);
set(h, 'EdgeColor', 'none');

colormap(flipud(cold)); 
cbr=colorbar;
cbr.TickLabelInterpreter='latex';
cbr.Ticks=[];
cbr.Title.Interpreter='latex';
cbr.Title.FontSize=24;
colorTitleHandle = get(cbr,'Title');
titleString = '$I(\textbf{r})$';
set(colorTitleHandle ,'String',titleString);
set(cbr,'position',[.915 .38 .02 .34])

set(h,'EdgeAlpha',0,'FaceAlpha',0.9);
caxis([0,1])
y1=[0,0];
x1=[xmin,xmax];
z1=[zmin,zmin];
plot3(z1+L,y1,x1,'Color',[0.7,0.7,0.7])
y2=[0,0];
x2=[xmin,xmax];
z2=[zmax,zmax];
plot3(z2+L,y2,x2,'Color',[0.7,0.7,0.7])
y3=[0,0];
x3=[xmin,xmin];
z3=[zmin,zmax];
plot3(z3+L,y3,x3,'Color',[0.7,0.7,0.7])
y4=[0,0];
x4=[xmax,xmax];
z4=[zmin,zmax];
plot3(z4+L,y4,x4,'Color',[0.7,0.7,0.7])
axis([zmin+L,zmax+L,xmin,xmax,xmin,xmax])
box on

rotate(h,[1,0,0],90)
view([-33,20])

set(gca,'FontSize',24)
ax = gca;

ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.ZLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';

xlabel('$z/\lambda_0$')
ylabel('$y/\lambda_0$')
zlabel('$x/\lambda_0$')

%% Compute temporal evolution

addpath('expmv-master')

tic

Deltac=0;

dt=1e-2;

g=sqrt(gammad*gammab)/(2*sqrt(2));

g1=sqrt(2)*g*gvec(1:na);
g2=sqrt(2)*g*gvec(na+1:end);

Totaltime=1.25*2*pi/g;

timescale=0:dt:Totaltime;


dim=2+2*na;

Psi_in=[1;zeros(dim-1,1)];

Hc=zeros(dim,dim);

Hc(1,1)=Hc(1,1)-Deltac;
Hc(2,2)=Hc(2,2)+Deltac;

Hg1=zeros(dim,dim);

for j=1:na
   Hg1(1,2+(j-1)+1)=conj(g1(j));
   Hg1(2+(j-1)+1,1)=g1(j); 
end

Hg2=zeros(dim,dim);

for j=1:na
   Hg2(2,2+na+(j-1)+1)=conj(g2(j));
   Hg2(2+na+(j-1)+1,2)=g2(j);
end

Hnh=zeros(dim,dim);

for j=1:2*na
    Hnh(2+(j-1)+1,2+(j-1)+1)=-1j*gamma/2.-1.*Deltalight;
    for k=j+1:2*na
       Jhop=1.5*gamma*lambda0*GreensTensor(lambda0,sqrt((xvec(k)-xvec(j))^2+(yvec(k)-yvec(j))^2),zvec(k)-zvec(j));
       Hnh(2+(k-1)+1,2+(j-1)+1)=-Jhop;
       Hnh(2+(j-1)+1,2+(k-1)+1)=-Jhop;
    end
end

H12=Hc+Hnh+Hg1+Hg2;

Probarray=zeros(size(timescale));
ProbQ1=zeros(size(timescale));
ProbQ2=zeros(size(timescale));

psi=Psi_in;

disp('evol')

for t=1:length(timescale)

    psi=expmv(-1j*dt,H12,psi);
    psiQ1=psi(1);
    psiQ2=psi(2);
    ProbQ1(t)=psiQ1'*psiQ1;
    ProbQ2(t)=psiQ2'*psiQ2;
    
    psirest=psi(3:end);
    Probarray(t)=psirest'*psirest;
end

%%

myBlue=[0.6,0.8,1];
myRed=[1,0.5,0.4];

figure
hold on
h=plot(timescale/(2*pi/g),Probarray,'--k');
h.LineWidth=4;
h=plot(timescale/(2*pi/g),ProbQ2);
h.LineWidth=4;
h.Color=myBlue;
h=plot(timescale/(2*pi/g),ProbQ1,'b');
h.LineWidth=4;
h.Color=myRed;

set(gca,'FontSize',30)
ax = gca;
grid on
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
xlabel('$\Omega t/2\pi$')
ylabel('Population')
ax.YScale='linear';
ax.XLim=[min(timescale/(2*pi/g)),max(timescale/(2*pi/g))];

%% FFT

kx=-pi/dperp+2*pi/(nperp*dperp)*(0:nperp-1);
ky=-pi/dperp+2*pi/(nperp*dperp)*(0:nperp-1);
kz=-pi/(2*L)+2*pi/(2*2*L)*(0:1);

ksq=zeros(2*nperp*nperp,1);
evenodd=zeros(2*nperp*nperp,1);
for j=1:2*na
    ueig=u(:,j);
    utens=reshape(ueig,[nperp,nperp,2]);
    ufft=fftn(utens)/sqrt(2*nperp*nperp);
    ufft=fftshift(ufft);
    ksq1=0;
    ksq2=0;
    for jx=1:nperp
        for jy=1:nperp
            ksq1=ksq1+sqrt(kx(jx)^2+ky(jy)^2)*abs(ufft(jx,jy,1))^2;
            ksq2=ksq2+sqrt(kx(jx)^2+ky(jy)^2)*abs(ufft(jx,jy,2))^2;
        end
    end
    if abs(ksq1-ksq2)/(ksq1+ksq2)<0.99
        error('ksq1 and ksq2 are both non-zero')
    end
    [ksq(j),evenodd(j)]=max([ksq1,ksq2]);
end

figure
hold on

for itera=1:2*na
    if evenodd(itera)==2
        h=scatter(ksq(itera)/(2*pi),-2*imag(lam(itera))/gamma,'w','filled','MarkerEdgeColor',0.*[1,1,1]);
        h.SizeData=200;
    else 
        h=scatter(ksq(itera)/(2*pi),-2*imag(lam(itera))/gamma,'filled','MarkerEdgeColor',0*[1,1,1],'MarkerFaceColor',0.45*[1,1,1]);
        h.SizeData=200;
    end
end
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'FontSize',36)
ax=gca;
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.YLabel.String='$\gamma_n/\gamma_e$';
ax.XLabel.String='$\overline{q}/k_0$';
ax.YTick=[1e-3,1e-2,1e-1,1e0];
ax.YLim=[1e-3,3e0];
ax.XLim=[2e-2,5e0];
ax.XTick=[1e-1,1e0];
box on
grid
ax.XMinorGrid='off';
ax.YMinorGrid='off';
grid off


%% Eigenvalues

figure
h=scatter(real(lam)/gamma,-2*imag(lam)/gamma);
set(gca,'yscale','log')
set(gca,'FontSize',36)
ax=gca;
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.XLabel.String='$\mathcal R(\epsilon_j)/\gamma_e$';
ax.YLabel.String='$-2 \mathcal I(\epsilon_j)/\gamma_e$';
ax.YTick=[1e-3,1e-2,1e-1,1e0];
ax.YLim=[1e-3,2e0];
ax.XLim=[-0.4,0.2];
h.SizeData=80;

%% Eigenmode

[sortlam,indexsort]=sort(imag(lam));
lamsort=lam(indexsort);
ueig=u(:,indexsort(end-3))-u(:,indexsort(end-4));

x_eig=xvec(1:na);
y_eig=yvec(1:na);
u_eig=ueig(1:na)*sqrt(2);


figure('Position', [10 10 150 150])
h=scatter(x_eig,y_eig,50,abs(u_eig),'filled');
h.MarkerEdgeColor=[1,1,1];
colormap(brewermap([],'Reds'))
caxis([min(abs(u_eig))-0.1*max(abs(u_eig)),1.25*max(abs(u_eig))]);
box on

set(gca,'FontSize',22)
ax = gca;
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.XTick=[];
ax.YTick=[];
ax.XLim=1.1*[min(xvec),max(xvec)];
ax.YLim=1.1*[min(xvec),max(xvec)];

figure('Position', [10 10 150 150])
h=scatter(x_eig,y_eig,50,angle(1j*u_eig),'filled');
h.MarkerEdgeColor=[1,1,1];
phasemap
box on

set(gca,'FontSize',22)
ax = gca;
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.XTick=[];
ax.YTick=[];
ax.XLim=1.1*[min(xvec),max(xvec)];
ax.YLim=1.1*[min(xvec),max(xvec)];

% figure
% colormap(flipud(hot)); 
% h=scatter(x_eig,y_eig,50,angle(u_eig),'filled');
% h.MarkerEdgeColor='k';
% col=colorbar;
% col.Ticks=[-pi,-pi/2,0,pi/2,pi];
% 
% set(gca,'FontSize',22)
% ax = gca;
% ax.FontName = 'LaTeX';
% ax.Title.Interpreter = 'LaTeX';
% ax.XLabel.Interpreter = 'LaTeX';
% ax.YLabel.Interpreter = 'LaTeX';
% ax.TickLabelInterpreter = 'LaTeX';
% xlabel('$x/\lambda_0$')
% ylabel('$y/\lambda_0$')


%% 

[sortlam,indexsort]=sort(imag(lam));
lamsort=lam(indexsort);
ueig=u(:,indexsort(end));

x_eig=xvec(1:na);
y_eig=yvec(1:na);
u_eig=ueig(1:na)*sqrt(2);


figure('Position', [10 10 150 150])
h=scatter(x_eig,y_eig,80,abs(u_eig),'filled','MarkerEdgeColor','w','LineWidth',1);
colormap(brewermap([],'Reds'))
caxis([min(abs(u_eig))-0.1*max(abs(u_eig)),1.25*max(abs(u_eig))]);
box off

set(gca,'FontSize',22)
ax = gca;
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.XTick=[];
ax.YTick=[];
ax.XLim=1.1*[min(xvec),max(xvec)];
ax.YLim=1.1*[min(xvec),max(xvec)];
box on

