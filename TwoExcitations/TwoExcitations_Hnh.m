clearvars
tic
addpath('expmv-master')

lambda0=1;
nperp=6;
dperp=0.5;
L=3;
gamma=1e1;
omega0=find_omega0opt(lambda0,nperp,dperp,L)

na=nperp^2;
zR=pi*omega0^2/lambda0;

xvec=setXpos(nperp,dperp);
yvec=setYpos(nperp,dperp);
zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
[gammad,gammab,Deltalight,overlaps]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

ratio=gammad/gammab

[gammad2,eigvec2,sortoverlaps,Deltalight2]=getShift2(gvec,lambda0,gamma,xvec,yvec,zvec);
toc
%%

ratio2=gammad2/(2*gammab);

Prob=zeros(1,2*na);

for j=1:2*na
    for k=j+1:2*na
        index= 2*na*(j-1)-j*(j-1)/2+(k-j-1);
        Prob(j)=Prob(j)+abs(eigvec2(index+1))^2;
        Prob(k)=Prob(k)+abs(eigvec2(index+1))^2;
    end
end
    
Probtrunc=Prob(1:na);
x_eig=xvec(1:na);
y_eig=yvec(1:na);

figure
colormap(flipud(hot)); 
h=scatter(x_eig,y_eig,50,Probtrunc,'filled');
h.MarkerEdgeColor='k';
col=colorbar;
col.Ticks=[];

set(gca,'FontSize',22)
ax = gca;
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
xlabel('$x/\lambda_0$')
ylabel('$y/\lambda_0$')

%%

Prob=zeros(1,2*na);

dxvec=0:dperp:(nperp-1)*dperp;
dyvec=0:dperp:(nperp-1)*dperp;
Pdxy=zeros(length(dxvec),length(dyvec));
rdxy=zeros(length(dxvec),length(dyvec));

for j=1:2*na
    for k=j+1:2*na
        index= 2*na*(j-1)-j*(j-1)/2+(k-j-1);
        dx=abs(xvec(j)-xvec(k));
        dy=abs(yvec(j)-yvec(k));
        [~,indx]=min(abs(dxvec-dx));
        [~,indy]=min(abs(dyvec-dy));
        Pdxy(indx,indy)=Pdxy(indx,indy)+abs(eigvec2(index+1))^2;
        rdxy(indx,indy)=sqrt(dxvec(indx)^2+dyvec(indy)^2);
    end
end

Pdxy(1,:)=Pdxy(1,:)*1;
Pdxy(:,1)=Pdxy(:,1)*1;
    
figure
h=contourf(Pdxy);
col=colorbar;
col.Ticks=[];


figure
h=contourf(rdxy);
col=colorbar;
col.Ticks=[];

rvec=reshape(rdxy,[1,nperp^2]);
Pvec=reshape(Pdxy,[1,nperp^2]);
figure
plot(rvec,Pvec,'o')

filter=0:dperp:2*nperp*dperp;
Pfilter=zeros(size(filter));

for iter=1:length(filter)-1
    Pfilter(iter)=sum(Pvec.*(rvec>=filter(iter)).*(rvec<filter(iter+1)))/(filter(iter+1)-filter(iter));
end

filters=filter+0.5*(filter(2)-filter(1));
figure
plot(filter,Pfilter)
hold on
plot(filter,filters.*exp(-filters.^2/(3*omega0^2)))
%%

Split=zeros(1,na);
NotSplit=zeros(1,na);

for j=1:na
    for k=j+1:na
        index= 2*na*(j-1)-j*(j-1)/2+(k-j-1);
        NotSplit(j)=NotSplit(j)+abs(eigvec2(index+1))^2;
        NotSplit(k)=NotSplit(k)+abs(eigvec2(index+1))^2;
    end
    for k=na+1:2*na
        index= 2*na*(j-1)-j*(j-1)/2+(k-j-1);
        Split(j)=Split(j)+abs(eigvec2(index+1))^2;
    end
end
x_eig=xvec(1:na);
y_eig=yvec(1:na);

figure
colormap(flipud(hot)); 
h=scatter(x_eig,y_eig,50,NotSplit,'filled');
h.MarkerEdgeColor='k';
col=colorbar;
col.Ticks=[];

set(gca,'FontSize',22)
ax = gca;
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
xlabel('$x/\lambda_0$')
ylabel('$y/\lambda_0$')

x_eig=xvec(1:na);
y_eig=yvec(1:na);

figure
colormap(flipud(hot)); 
h=scatter(x_eig,y_eig,50,Split,'filled');
h.MarkerEdgeColor='k';
col=colorbar;
col.Ticks=[];

set(gca,'FontSize',22)
ax = gca;
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
xlabel('$x/\lambda_0$')
ylabel('$y/\lambda_0$')

%%

x_eig=xvec(1:na);
y_eig=yvec(1:na);

figure
for j=1:na
    subplot(nperp,nperp,floor((na-j)/nperp)*nperp+mod(j-1,nperp)+1)
    condProb=zeros(1,na);
    for k=na+1:2*na
        index= 2*na*(j-1)-j*(j-1)/2+(k-j-1);
        condProb(k-na)=abs(eigvec2(index+1))^2;
    end
    h=scatter(x_eig,y_eig,50,(condProb/max(abs(eigvec2).^2)).^(0.5),'filled');
    caxis([0,1])    
    h.MarkerEdgeColor='k';

    ax = gca;
    ax.XTick=[];
    ax.YTick=[];
end
colormap(flipud(hot)); 

figure
for j=1:na
    subplot(nperp,nperp,floor((na-j)/nperp)*nperp+mod(j-1,nperp)+1)
    condProb=zeros(1,na);
    for k=1:j-1
        index= 2*na*(k-1)-k*(k-1)/2+(j-k-1);
        condProb(k)=abs(eigvec2(index+1))^2;
    end
    for k=j+1:na
        index= 2*na*(j-1)-j*(j-1)/2+(k-j-1);
        condProb(k)=abs(eigvec2(index+1))^2;
    end
    h=scatter(x_eig,y_eig,50,(condProb/max(abs(eigvec2).^2)).^(0.5),'filled');
    caxis([0,1])
    h.MarkerEdgeColor='k';

    ax = gca;
    ax.XTick=[];
    ax.YTick=[];
end
colormap(flipud(hot)); 

%%

x_eig=xvec(1:na);
y_eig=yvec(1:na);

figure
for j=1:na
    subplot(nperp,nperp,floor((na-j)/nperp)*nperp+mod(j-1,nperp)+1)
    condProb=zeros(1,na);
    for k=na+1:2*na
        index= 2*na*(j-1)-j*(j-1)/2+(k-j-1);
        condProb(k-na)=abs(eigvec2(index+1))^2;
    end
    h=scatter(x_eig,y_eig,50,(condProb),'filled');
    caxis([0,max((condProb))])    
    h.MarkerEdgeColor='k';

    ax = gca;
    ax.XTick=[];
    ax.YTick=[];
end
colormap(flipud(hot)); 

figure
for j=1:na
    subplot(nperp,nperp,floor((na-j)/nperp)*nperp+mod(j-1,nperp)+1)
    condProb=zeros(1,na);
    for k=1:j-1
        index= 2*na*(k-1)-k*(k-1)/2+(j-k-1);
        condProb(k)=abs(eigvec2(index+1))^2;
    end
    for k=j+1:na
        index= 2*na*(j-1)-j*(j-1)/2+(k-j-1);
        condProb(k)=abs(eigvec2(index+1))^2;
    end
    h=scatter(x_eig,y_eig,50,(condProb),'filled');
    caxis([0,max(condProb)])
    h.MarkerEdgeColor='k';

    ax = gca;
    ax.XTick=[];
    ax.YTick=[];
end
colormap(flipud(hot)); 
