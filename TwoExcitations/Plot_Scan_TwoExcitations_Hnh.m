
load('2exc_scan.mat')
load('gvec_mean.mat')

gamma=10;

nperpvec=[4,5,6,7,8];
Nnperp=length(nperpvec);

figure

myBlue=[0.6,0.8,1];
plot([0,1],[0,1],'w');
hold on
for iternperp=1:Nnperp
    xi=xioptmat(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    ratio=ratiomat(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    ratio2=ratiomat2(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    h=plot(xi,ratio2./ratio,'.-');
    h.MarkerSize=25;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    hold on
end

set(gca,'FontSize',22)
ax = gca;
ax.XLim = [1e-1, 1e1];
ax.XTick=[1e-1,1e0,1e1];
ax.XScale='log';
ax.YScale='log';

ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
grid on
xlabel('$L\lambda_0/L_\perp^2$')
ylabel('$\gamma_d^{(2)}/(2\gamma_d)$')


%%
figure

myBlue=[0.6,0.8,1];
plot([0,1],[0,1],'w');
hold on
for iternperp=1:Nnperp
    gb=gammabvec(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>.1);
    xi=(0.5*nperpvec(iternperp))^2.*xioptmat(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>.1);
    gd=gb.*ratiomat(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>.1);
    gd2=(2*gb).*ratiomat2(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>.1);
    h=plot(xi,gd2./(2*gd),'.-');
    h.MarkerSize=25;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    hold on
    
    
    vmean=gvec_mean(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>.1);
    gd_th=2*gd.*(1-2*vmean)+2*vmean*gamma;
    h=plot(xi,gd_th./(2*gd),'--');
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
end

set(gca,'FontSize',22)
ax = gca;
ax.XScale='linear';
ax.YScale='log';
ax.XLim=[1e0,1e1];

ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
grid on
xlabel('$L/\lambda_0$')
ylabel('$\gamma_d^{(2)}/(2\gamma_d)$')
