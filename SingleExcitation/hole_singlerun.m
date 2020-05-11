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

%% 1 hole

iterhole=round(na+nperp*round(nperp/2)+nperp/2);
iterhole=105;

xvec=setXpos(nperp,dperp);
yvec=setYpos(nperp,dperp);
zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);

Hnh=zeros(length(xvec),length(xvec));
        
for j=1:length(xvec)
    Hnh(j,j)=Hnh(j,j)-(1-(j==iterhole))*1j*gamma/2.;
    for k=j+1:length(xvec)
       Jhop=(1-(k==iterhole))*(1-(j==iterhole))*1.5*gamma*lambda0*GreensTensor(lambda0,sqrt((xvec(k)-xvec(j))^2+(yvec(k)-yvec(j))^2),zvec(k)-zvec(j));
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

gammad;
gammab;
ratiovec=gammad/gammab

%% 2holes

iterhole2=round(nperp*round(nperp/2)+nperp/2);

xvec(iterhole2)=xvec(iterhole2)+(nperp*dperp)*10000;

gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
[gammad,gammab,Deltalight,overlaps,Deltabright,ueig]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

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
