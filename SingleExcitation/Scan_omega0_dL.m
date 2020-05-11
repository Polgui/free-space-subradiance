clearvars

tic

lambda0=1;

omega0min=0.5*lambda0;
omega0max=5*lambda0;
Nomega0=100;
omega0vec=linspace(omega0min,omega0max,Nomega0);

dLmin=-0.25;
dLmax=0.25;
NdL=201;
dLvec=linspace(dLmin,dLmax,NdL);

Lav=15;

omega0opt=0;
Lopt=0;
        
nperp=12;
dperp=0.8;
gamma=1e1;
ratiovec=1;

        
for iteromega0=1:Nomega0
    iteromega0/Nomega0
    
    omega0=omega0vec(iteromega0);
    zR=pi*omega0^2/lambda0;

    xvec=setXpos(nperp,dperp);
    yvec=setYpos(nperp,dperp);
    zvec=setZpos(nperp,Lav,lambda0,zR,xvec,yvec);

    gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
    [gammad,gammab,Deltalight,overlaps]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

    if gammad/gammab<ratiovec
        ratiovec=gammad/gammab;
        omega0opt=omega0;
    end

end

omega0=omega0opt;
zR=pi*omega0^2/lambda0;
    
ratiovec=ones(1,length(dLvec));
Shiftvec=ones(1,length(dLvec));
ueigvat=ones(2*nperp^2,length(dLvec));
ueigvat2=ones(2*nperp^2,length(dLvec));

for iterL=1:NdL
    iterL/NdL

    xvec=setXpos(nperp,dperp);
    yvec=setYpos(nperp,dperp);
    zvec=setZpos(nperp,Lav,lambda0,zR,xvec,yvec,dLvec(iterL));

    gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
    [gammad,gammab,Deltalight,overlaps,Deltabright,ueig,ueig2]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

    ratiovec(iterL)=gammad/gammab;
    Shiftvec(iterL)=(Deltalight-Deltabright)/gamma;
    
    ueigvat(:,iterL)=ueig;
    ueigvat2(:,iterL)=ueig2;
end
toc

%% Plot

Lvec=2*(Lav+dLvec);
Shiftvec_new=Shiftvec;
ratiovec_new=ratiovec;

Shiftvec_new(Lvec<30-0.25)=-Shiftvec(Lvec<30-0.25);
Shiftvec_new(Lvec>30+0.25)=-Shiftvec(Lvec>30+0.25);
ratiovec_new(Lvec<30-0.25)=1./ratiovec(Lvec<30-0.25);
ratiovec_new(Lvec>30+0.25)=1./ratiovec(Lvec>30+0.25);

figure('pos',[10 10 500 200])
yyaxis left
h=plot(Lvec,Shiftvec_new,'--', 'linewidth',2);
set(gca,'FontSize',30)
ax = gca;
ax.YLabel.Interpreter = 'LaTeX';
ylabel('$(\Delta_a-\Delta_s)/\gamma_e$')
yyaxis right 
h=semilogy(Lvec,ratiovec_new,'linewidth',2);
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
ax.YTick=[1e-4,1e-2,1e0,1e2,1e4];
xlabel('$L/\lambda_0$')
ylabel('$\gamma_a/\gamma_s$')

%%
% line1.LineStyle='--';
% line1.LineWidth=4;
% line2.Color='r';
% line2.LineWidth=4;
% AX(1).Position=[0.1650    0.220    0.6500    0.7150];
% AX(2).YScale='linear';
% AX(2).YAxis.Color='r';
% AX(1).LineWidth=1;
% AX(2).LineWidth=1;
% AX(1).TickLength=[0.02,0.05];
% AX(2).TickLength=[0.02,0.05];
% AX(1).XLabel.String='$L/\lambda_0$';
% AX(1).YLabel.String='$(\Delta_a-\Delta_s)/\gamma_e$';
% AX(2).YLabel.String='log${}_{10}(\gamma_a/\gamma_s)$';
% 
% AX(1).set('FontSize',30);
% AX(1).YLabel.Position=[29.3892    -0.1500   -1.0000];
% AX(1).FontName = 'LaTeX';
% AX(1).Title.Interpreter = 'LaTeX';
% AX(1).XLabel.Interpreter = 'LaTeX';
% AX(1).YLabel.Interpreter = 'LaTeX';
% AX(1).TickLabelInterpreter = 'LaTeX';
% 
% 
% %AX(2).YLabel.Position=[30.5542 0.0316 -1];
% 
% 
% AX(2).set('FontSize',30);
% AX(2).FontName = 'LaTeX';
% AX(2).Title.Interpreter = 'LaTeX';
% AX(2).XLabel.Interpreter = 'LaTeX';
% AX(2).YLabel.Interpreter = 'LaTeX';
% AX(2).TickLabelInterpreter = 'LaTeX';
% 
% %AX(2).YTick=[-4,-2,0];
% %AX(2).YLim=[-4,0];
% 
% figure('pos',[10 10 500 200])
% [AX,line1,line2]=plotyy(Lvec,Shiftvec_new,Lvec,log10(ratiovec_new));
% line1.LineStyle='--';
% line1.LineWidth=4;
% line2.Color='r';
% line2.LineWidth=4;
% AX(1).Position=[0.1650    0.220    0.6500    0.7150];
% AX(2).YScale='linear';
% AX(2).YAxis.Color='r';
% AX(1).LineWidth=1;
% AX(2).LineWidth=1;
% AX(1).TickLength=[0.02,0.05];
% AX(2).TickLength=[0.02,0.05];
% AX(1).XLabel.String='$L/\lambda_0$';
% AX(1).YLabel.String='$(\Delta_a-\Delta_s)/\gamma_e$';
% AX(2).YLabel.String='log${}_{10}(\gamma_a/\gamma_s)$';
% 
% AX(1).set('FontSize',30);
% AX(1).YLabel.Position=[29.3892    -0.1500   -1.0000];
% AX(1).FontName = 'LaTeX';
% AX(1).Title.Interpreter = 'LaTeX';
% AX(1).XLabel.Interpreter = 'LaTeX';
% AX(1).YLabel.Interpreter = 'LaTeX';
% AX(1).TickLabelInterpreter = 'LaTeX';
% 
% 
% %AX(2).YLabel.Position=[30.5542 0.0316 -1];
% 
% 
% AX(2).set('FontSize',30);
% AX(2).FontName = 'LaTeX';
% AX(2).Title.Interpreter = 'LaTeX';
% AX(2).XLabel.Interpreter = 'LaTeX';
% AX(2).YLabel.Interpreter = 'LaTeX';
% AX(2).TickLabelInterpreter = 'LaTeX';