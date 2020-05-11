clearvars

load('Compress.mat')

nperpvec=[4,8,12,16,20];
Nperp=5;

lambda0=1;

dperpvec=linspace(0.1,1.,24);
Ndperp=length(dperpvec);
        
gamma=1e1;

Delta_vec=zeros(Nperp,Ndperp);
Lvec=zeros(Nperp,Ndperp);
ratiovec=zeros(Nperp,Ndperp);


mycluster=parcluster();
mypool=parpool(mycluster); 

for iterperp=1:Nperp
    nperp=nperpvec(iterperp);
    parfor iterdperp=1:Ndperp

        dperp=dperpvec(iterdperp);

        L=1.25;

        omega0=find_omega0opt(lambda0,nperp,dperp,L);
        zR=pi*omega0^2/lambda0;

        xvec=setXpos(nperp,dperp);
        yvec=setYpos(nperp,dperp);
        zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);
        Lvec(iterperp,iterdperp)=min(zvec);

        gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
        [gammad,gammab,Deltalight,overlaps]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

        Delta_vec(iterperp,iterdperp)=Deltalight/gamma;
        ratiovec(iterperp,iterdperp)=gammad/gammab;

    end
end

delete(gcp)

save('shortL.mat','ratiovec','Delta_vec','Lvec','dperpvec','nperpvec')
%%
load('Compress')
load('shortL')

Nperp=5;
Ndperp=length(dperpvec);
        
gamma=1e1;
figure
subplot(2,1,1)

deltaperpvec=0.1:0.1:1.0;
nperpvec=[4,8,12,16,20];
Nnperp=length(nperpvec);

myBlue=[0.6,0.8,1];

for iternperp=1:length(nperpvec)
    ratiomat=zeros(1,length(deltaperpvec));
    for iterdperp=1:1:10
        ratiovec=ratiovec_cell{iterdperp};
        xioptvec=xiopt_cell{iterdperp};

        ratio=ratiovec(iternperp,:);
        xiopt=xioptvec(iternperp,:);
        ratiocut=ratio(xiopt.*(0.1*iterdperp*nperpvec(iternperp))^2>1);
        ratiomat(iterdperp)=min(ratiocut);

    end
    h=semilogy(deltaperpvec,ratiomat*nperpvec(iternperp)^4,'.-');
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    h.LineWidth=6.0;
    h.MarkerSize=45;
    hold on
end

set(gca,'FontSize',30)
ax = gca;
%ax.Position=[0.1447,3*0.1501,0.7603,0.7749/2];
ax.FontName = 'LaTeX';
grid on
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.ZLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.YLim=[1e-2,1e4];
ax.YTick=[1e-2,1e0,1e2,1e4];
ax.XTick=[];
ylabel('$n_\perp^4 \gamma_d/\gamma_b$')


subplot(2,1,2)

for iterperp=1:Nperp
    h=plot(dperpvec(1:end-1),Delta_vec(iterperp,1:end-1),'.-');    
    hold on
    h.MarkerSize=45;
    h.LineWidth=6;
    h.Color=myBlue*(5-iterperp+1)/5;
end

set(gca,'FontSize',30)
ax = gca;

ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.YLim=[-3,1];
grid on
xlabel('$\delta_\perp/\lambda_0$')
ylabel('$\Delta_d/\gamma$')
