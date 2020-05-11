clearvars

tic

nperp=12;
dperp=0.8;
L=15;
2*L/(nperp*dperp)^2
gamma=1e1;

na=nperp^2;


lambda0=1;
omega0min=0.5*lambda0;
omega0max=5*lambda0;
Nomega0=100;
omega0vec=linspace(omega0min,omega0max,Nomega0);

ratiovec=zeros(1,Nomega0);

minzvec=zeros(1,Nomega0);
maxzvec=zeros(1,Nomega0);
allzvec=zeros(2*na,Nomega0);

n0vec=zeros(1,Nomega0);

for iteromega0=1:Nomega0

    omega0=omega0vec(iteromega0);
    zR=pi*omega0^2/lambda0;


    xvec=setXpos(nperp,dperp);
    yvec=setYpos(nperp,dperp);
    [zvec,n0]=setZpos(nperp,L,lambda0,zR,xvec,yvec);

    gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
    [gammad,gammab,Deltalight,overlaps]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

    ratiovec(iteromega0)=gammad/gammab;
    
    minzvec(iteromega0)=min(abs(zvec));
    maxzvec(iteromega0)=max(abs(zvec));
    
    allzvec(:,iteromega0)=abs(zvec);
    
    n0vec(iteromega0)=n0;

end
toc

%% Plot


figure
[AX,line1,line2]=plotyy([omega0vec',omega0vec'],[minzvec',maxzvec'],omega0vec,ratiovec);
line1(1).Color='k';
line1(1).LineStyle='--';
line1(1).LineWidth=2;
line1(2).Color='k';
line1(2).LineStyle='--';
line1(2).LineWidth=2;
line2.Color='r';
line2.LineWidth=2;
AX(1).Position=[0.1050    0.1500    0.7500    0.8150];
AX(2).YScale='log';
AX(2).YAxis.Color='r';
AX(1).LineWidth=1;
AX(2).LineWidth=1;
AX(1).TickLength=[0.02,0.05];
AX(2).TickLength=[0.02,0.05];
AX(1).XLabel.String='$w_0/\lambda_0$';
AX(1).YLabel.String='$z/\lambda_0$';
AX(2).YLabel.String='$\gamma_d/\gamma_b$';

AX(1).set('FontSize',22);
AX(1).FontName = 'LaTeX';
AX(1).Title.Interpreter = 'LaTeX';
AX(1).XLabel.Interpreter = 'LaTeX';
AX(1).YLabel.Interpreter = 'LaTeX';
AX(1).TickLabelInterpreter = 'LaTeX';


AX(2).set('FontSize',22);
AX(2).FontName = 'LaTeX';
AX(2).Title.Interpreter = 'LaTeX';
AX(2).XLabel.Interpreter = 'LaTeX';
AX(2).YLabel.Interpreter = 'LaTeX';
AX(2).TickLabelInterpreter = 'LaTeX';

AX(2).YTick=[1e-3,1e-2,1e-1];
AX(2).YLim=[1e-3,1e-1];

% R= maxzvec.*(1+(pi*omega0vec.^2./maxzvec).^2);
% thetamax=atan(nperp*dperp./(sqrt(2)*R));
% z=cos(thetamax).*R;
% dz=maxzvec-z;

%figure
%plot(omega0vec,R)
%figure
%plot(omega0vec,thetamax)
%figure
%plot(omega0vec,z)


%figure
%plot(omega0vec,n0vec)


%figure
%plot(omega0vec,allzvec)


% omegamaxth=-gamma*sin(phivec);
% figure
% h=semilogx(phivec/pi,omegamaxth/2,'--r');
% h.LineWidth=2;
% hold on
% for itergamma=length(gammapvec):-1:1
%     h=semilogx(phivec(imag(deltadark(itergamma,:))==0)/pi,omegamax(itergamma,imag(deltadark(itergamma,:))==0)/2,'.-');
%     h.MarkerSize=25;
%     h.LineWidth=1;
%     h.Color=myBlue*itergamma/length(gammapvec);
%     hold on
%     
% end
% 
% set(gca,'FontSize',22)
% ax = gca;
% ax.YLim=[-0.2,0.];
% ax.YTick=(-0.2:0.05:0);
% 
% ax.FontName = 'LaTeX';
% ax.Title.Interpreter = 'LaTeX';
% ax.XLabel.Interpreter = 'LaTeX';
% ax.YLabel.Interpreter = 'LaTeX';
% ax.TickLabelInterpreter = 'LaTeX';
% 
% xlabel('${\phi}/{\pi}$')
% ylabel('$\omega_{p}/\gamma_{1D}$')
