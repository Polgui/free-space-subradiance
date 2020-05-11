%% gammad/gammab
clearvars


lambda0=1;


nperp=6;
dperp=0.5;
L=2;
gamma=1e1;

na=nperp^2;
dperp*nperp
omega0=find_omega0opt(lambda0,nperp,dperp,L)
zR=pi*omega0^2/lambda0;

xvec=setXpos(nperp,dperp);
yvec=setYpos(nperp,dperp);
zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

%%

[gvec,GaussianBeam]=Hopping(xvec,yvec,zvec,lambda0,omega0);
[gammad,gammab,Deltalight,overlaps,Deltabright,ueig,ueig2,lam,u]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

%% 

Deltavec=linspace(-0.4*gamma,0.4*gamma,501);
Rplus=zeros(1,length(Deltavec));
Rminus=zeros(1,length(Deltavec));

for iterDelta=1:length(Deltavec)
    Delta=Deltavec(iterDelta);

    Hnh=zeros(length(xvec),length(xvec));

    for j=1:length(xvec)
        Hnh(j,j)=Hnh(j,j)-Deltalight+Delta-1j*gamma/2.;
        for k=j+1:length(xvec)
           Jhop=1.5*gamma*lambda0*GreensTensor(lambda0,sqrt((xvec(k)-xvec(j))^2+(yvec(k)-yvec(j))^2),zvec(k)-zvec(j));
           Hnh(k,j)=Hnh(k,j)-Jhop;
           Hnh(j,k)=Hnh(j,k)-Jhop;
        end
    end

    Hinv=inv(Hnh);

    Rplus(iterDelta)=abs(GaussianBeam*Hinv*transpose(GaussianBeam))^2*9*pi^2*gamma^2/(4*(2*pi)^4);
    
    Hnh=zeros(length(xvec),length(xvec));

    for j=1:nperp*nperp
        Hnh(j,j)=Hnh(j,j)-Deltalight+Delta-1j*gamma/2.;
        for k=j+1:length(xvec)
           Jhop=1.5*gamma*lambda0*GreensTensor(lambda0,sqrt((xvec(k)-xvec(j))^2+(yvec(k)-yvec(j))^2),zvec(k)-zvec(j));
           Hnh(k,j)=Hnh(k,j)-Jhop;
           Hnh(j,k)=Hnh(j,k)-Jhop;
        end
    end
    for j=nperp*nperp+1:length(xvec)
        Hnh(j,j)=Hnh(j,j)-Deltalight-Delta-1j*gamma/2.;
        for k=j+1:length(xvec)
           Jhop=1.5*gamma*lambda0*GreensTensor(lambda0,sqrt((xvec(k)-xvec(j))^2+(yvec(k)-yvec(j))^2),zvec(k)-zvec(j));
           Hnh(k,j)=Hnh(k,j)-Jhop;
           Hnh(j,k)=Hnh(j,k)-Jhop;
        end
    end
    
    Hinv=inv(Hnh);
    Rminus(iterDelta)=abs(GaussianBeam*Hinv*transpose(GaussianBeam))^2*9*pi^2*gamma^2/(4*(2*pi)^4);
end

%%

Hnh=zeros(length(xvec),length(xvec));

for j=1:length(xvec)
    Hnh(j,j)=Hnh(j,j)-Deltalight-1j*gamma/2.;
    for k=j+1:length(xvec)
       Jhop=1.5*gamma*lambda0*GreensTensor(lambda0,sqrt((xvec(k)-xvec(j))^2+(yvec(k)-yvec(j))^2),zvec(k)-zvec(j));
       Hnh(k,j)=Hnh(k,j)-Jhop;
       Hnh(j,k)=Hnh(j,k)-Jhop;
    end
end

Hinv=inv(Hnh);
[U,V]=eig(Hinv);
GaussianBeamnorm=GaussianBeam/norm(GaussianBeam);
overlaps=GaussianBeamnorm*U;

[sortoverlaps,ind]=sort(abs(overlaps).^2);
E=Hinv(ind(end))
figure
plot(sortoverlaps,'x')

%%

myBlue=[0.6,0.8,1];
myRed=[1,0.5,0.4];

Gamma=3*pi*gamma/(2*pi*dperp)^2;
figure
h=plot(Deltavec/Gamma,Rplus);
h.LineWidth=4;
h.Color=myBlue;
hold on
h=plot(Deltavec/Gamma,Rminus);
h.Color=myRed;
h.LineWidth=4;
h=plot(Deltavec/Gamma,1./((1+gammad/(gammab-gammad))^2+4*(Deltavec).^2/(gammab-gammad)^2),'--k');
h.LineWidth=4;
h=plot(Deltavec/Gamma,(1./(1+gammad/(gammab-gammad)+4*Deltavec.^2/(gammad*(gammab-gammad)))).^2,'--k');
h.LineWidth=4;

set(gca,'FontSize',30)
ax = gca;
grid on
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
xlabel('$\Delta/\Gamma$')
ylabel('$R(\Delta)$')
ax.YScale='linear';
ax.XLim=[-1,1];
ax.XTick=-1:0.5:1;