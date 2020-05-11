%% Compress

clearvars

deltaperpvec=0.1:0.1:1.0;
nperpvec=[4,8,12,16,20];

omega0opt_cell=cell(1,10);
ratiovec_cell=cell(1,10);
xiopt_cell=cell(1,10);

load('27-Jun-2018 00:39:20.mat')
omega0opt_cell{1}=omega0opt;
ratiovec_cell{1}=ratiovec;
xiopt_cell{1}=xiopt;

load('27-Jun-2018 00:18:11.mat')
omega0opt_cell{2}=omega0opt;
ratiovec_cell{2}=ratiovec;
xiopt_cell{2}=xiopt;

load('27-Jun-2018 00:31:23.mat')
omega0opt_cell{3}=omega0opt;
ratiovec_cell{3}=ratiovec;
xiopt_cell{3}=xiopt;

load('27-Jun-2018 07:09:38.mat')
omega0opt_cell{4}=omega0opt;
ratiovec_cell{4}=ratiovec;
xiopt_cell{4}=xiopt;

load('27-Jun-2018 06:17:04.mat')
omega0opt_cell{5}=omega0opt;
ratiovec_cell{5}=ratiovec;
xiopt_cell{5}=xiopt;

load('27-Jun-2018 06:28:30.mat')
omega0opt_cell{6}=omega0opt;
ratiovec_cell{6}=ratiovec;
xiopt_cell{6}=xiopt;

load('27-Jun-2018 19:47:31.mat')
omega0opt_cell{7}=omega0opt;
ratiovec_cell{7}=ratiovec;
xiopt_cell{7}=xiopt;

load('27-Jun-2018 18:47:19.mat')
omega0opt_cell{8}=omega0opt;
ratiovec_cell{8}=ratiovec;
xiopt_cell{8}=xiopt;

load('27-Jun-2018 19:24:07.mat')
omega0opt_cell{9}=omega0opt;
ratiovec_cell{9}=ratiovec;
xiopt_cell{9}=xiopt;

load('27-Jun-2018 19:06:53.mat')
omega0opt_cell{10}=omega0opt;
ratiovec_cell{10}=ratiovec;
xiopt_cell{10}=xiopt;

save('Compress.mat','omega0opt_cell','ratiovec_cell','xiopt_cell')


   
%% Plots
close all

myBlue=[0.6,0.8,1];



figure
    
for iterdperp=1:1:10
    ratiovec=ratiovec_cell{iterdperp};
    xiopt=xiopt_cell{iterdperp};

    Nnperp=length(nperpvec);
    
    for iternperp=1:Nnperp
        h=plot3(xiopt(iternperp,xiopt(iternperp,:)>1e-2),deltaperpvec(iterdperp)*ones(size(xiopt(iternperp,xiopt(iternperp,:)>1e-2))),ratiovec(iternperp,xiopt(iternperp,:)>1e-2),'.-');
        h.MarkerSize=25;
        h.LineWidth=2;
        h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
        hold on
    end
end
set(gca,'FontSize',22)
ax = gca;
ax.XLim = [1e-2, 1e1];
ax.XTick=[1e-2,1e-1,1e0,1e1];
ax.XScale='log';
ax.ZScale='log';

ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.ZLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
grid on
xlabel('$L\lambda_0/L_\perp^2$')
ylabel('$\delta_\perp/\lambda_0$')
zlabel('$\gamma_d/\gamma_b$')

%% 
load('2exc_scan.mat')
load('gvec_mean.mat')

gamma=1e1;
figure
plot([0,1],[0,1],'w')
hold on
nperpvec=[4,5,6,7,8];
Nnperp=length(nperpvec);

myRed=[1,0.4,0.4];

for iternperp=1:2:Nnperp
    xi=xioptmat(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    ratio=ratiomat2(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    h=plot(xi(1:2:end),ratio(1:2:end),'.-','HandleVisibility','off');
    h.MarkerSize=30;
    h.LineWidth=2;
    h.Color=myRed*(Nnperp-iternperp+5)/(Nnperp+4);
    hold on
end

set(gca,'FontSize',30)
ax = gca;
ax.XLim = [1e-2, 1e1];
ax.XTick=[1e-2,1e-1,1e0,1e1];
ax.YTick=[1e-4,1e-2,1e0];
ax.XScale='log';
ax.YScale='log';

ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
grid
ax.XMinorGrid='off';
ax.YMinorGrid='off';
xlabel('$L\lambda_0/L_\perp^2$')
ylabel('$\gamma_d/\gamma_b$')

%mylegend=["2 exc., $N_\perp=4$";"2 exc., $N_\perp=6$";"2 exc., $N_\perp=8$"];
% hl=legend(mylegend);
% hl.Interpreter='LaTeX';
       
load('Compress.mat')
nperpvec=[4,8,12,16,20];

myBlue=[0.6,0.8,1];
iterdperp=5;
ratiovec=ratiovec_cell{iterdperp};
xiopt=xiopt_cell{iterdperp};

Nnperp=length(nperpvec);

for iternperp=1:1:Nnperp
    xi=xiopt(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    ratio=ratiovec(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    
    h=plot(xi,ratio,'.-');
    h.MarkerSize=30;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    hold on
end


mylegend=["1 exc."; "$N_\perp=4$";"$N_\perp=8$";"$N_\perp=12$";"$N_\perp=16$";"$N_\perp=20$"];
 hl=legend(mylegend,'FontSize',22);
 hl.Interpreter='LaTeX';
 hl.Position=[0.7703 0.4463 0.1172 0.3549];

myRed=[1,0.4,0.4];

for iternperp=1:2:Nnperp
    xi=xiopt_mean(iternperp,0.5*xiopt_mean(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    
    ratio=ratiomat(iternperp,0.5*xiopt_mean(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    vmean=gvec_mean(iternperp,0.5*xiopt_mean(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    gb=gammabvec(iternperp,0.5*xiopt_mean(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    gd=ratio.*gb;
    gd_th=gd.*(1-2*vmean)+vmean*gamma;

    h=plot(xi,(gd_th)./(gb),'--','HandleVisibility','off');
    h.MarkerSize=30;
    h.LineWidth=2;
    h.Color=myRed*(Nnperp-iternperp+8)/(Nnperp+7);
end
 
 
 xi=xiopt(1,0.5*xiopt(1,:).*(0.1*iterdperp*nperpvec(1))^2>1);
    ratio=ratiovec(1,0.5*xiopt(1,:).*(0.1*iterdperp*nperpvec(1))^2>1);
    
    h=plot(xi,ratio,'.-','HandleVisibility','off');
    h.MarkerSize=30;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-1+1)/Nnperp;
   
    
%% 
load('2exc_scan.mat')
load('gvec_mean.mat')

gamma=1e1;
figure
plot([0,1],[0,1],'w')
hold on
nperpvec=[4,5,6,7,8];
Nnperp=length(nperpvec);

myRed=[1,0.4,0.4];

for iternperp=1:2:Nnperp
    xi=xioptmat(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    ratio=ratiomat2(iternperp,0.5*xioptmat(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    h=plot(xi(1:2:end),ratio(1:2:end),'.-');
    h.MarkerSize=30;
    h.LineWidth=2;
    h.Color=myRed*(Nnperp-iternperp+5)/(Nnperp+4);
    hold on
end

set(gca,'FontSize',30)
ax = gca;
ax.XLim = [1e-2, 1e1];
ax.XTick=[1e-2,1e-1,1e0,1e1];
ax.YTick=[1e-4,1e-2,1e0];
ax.XScale='log';
ax.YScale='log';

ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
grid
ax.XMinorGrid='off';
ax.YMinorGrid='off';
xlabel('$L\lambda_0/L_\perp^2$')
ylabel('$\gamma_d/\gamma_b$')

%mylegend=["2 exc., $N_\perp=4$";"2 exc., $N_\perp=6$";"2 exc., $N_\perp=8$"];
% hl=legend(mylegend);
% hl.Interpreter='LaTeX';
       
load('Compress.mat')
nperpvec=[4,8,12,16,20];

myBlue=[0.6,0.8,1];
iterdperp=5;
ratiovec=ratiovec_cell{iterdperp};
xiopt=xiopt_cell{iterdperp};

Nnperp=length(nperpvec);

for iternperp=1:1:Nnperp
    xi=xiopt(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    ratio=ratiovec(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    
    h=plot(xi,ratio,'.-','HandleVisibility','off');
    h.MarkerSize=30;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    hold on
end


mylegend=["2 exc."; "$N_\perp=4$";"$N_\perp=6$";"$N_\perp=8$"];
 hl=legend(mylegend,'FontSize',22);
 hl.Interpreter='LaTeX';
 hl.Position=[0.7703 0.4463 0.1172 0.3549];

myRed=[1,0.4,0.4];

for iternperp=1:2:Nnperp
    xi=xiopt_mean(iternperp,0.5*xiopt_mean(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    
    ratio=ratiomat(iternperp,0.5*xiopt_mean(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    vmean=gvec_mean(iternperp,0.5*xiopt_mean(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    gb=gammabvec(iternperp,0.5*xiopt_mean(iternperp,:).*(0.5*nperpvec(iternperp))^2>1);
    gd=ratio.*gb;
    gd_th=gd.*(1-2*vmean)+vmean*gamma;

    h=plot(xi,(gd_th)./(gb),'--','HandleVisibility','off');
    h.MarkerSize=25;
    h.LineWidth=2;
    h.Color=myRed*(Nnperp-iternperp+8)/(Nnperp+7);
end
 
 
 xi=xiopt(1,0.5*xiopt(1,:).*(0.1*iterdperp*nperpvec(1))^2>1);
    ratio=ratiovec(1,0.5*xiopt(1,:).*(0.1*iterdperp*nperpvec(1))^2>1);
    
    h=plot(xi,ratio,'.-','HandleVisibility','off');
    h.MarkerSize=30;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-1+1)/Nnperp;
    
%%

% myBlue=[0.6,0.8,1];
% 
% figure
% ratiomat=zeros(length(nperpvec),length(deltaperpvec));
% for iterdperp=1:1:10
%     ratiovec=ratiovec_cell{iterdperp};
%     xioptvec=xiopt_cell{iterdperp};
%     for iternperp=1:length(nperpvec)
%         ratio=ratiovec(iternperp,:);
%         xiopt=xioptvec(iternperp,:);
%         ratiocut=ratio(xiopt.*(0.1*iterdperp*nperpvec(iternperp))^2>1);
%         ratiomat(iternperp,iterdperp)=min(ratiocut);
%     end
%     h=loglog(nperpvec,ratiomat(:,iterdperp));
%     h.Color=myBlue*(10-iterdperp+1)/10;
%     h.LineWidth=2.0;
%     hold on
% end
% h=loglog(nperpvec,2./nperpvec.^4,'--r');
% h.LineWidth=2.0;
% 
% set(gca,'FontSize',22)
% ax = gca;
% ax.Position=[0.1447,3*0.1501,0.7603,0.7749/2];
% ax.FontName = 'LaTeX';
% ax.XTickMode='manual';
% ax.XTick=[5,10,15,20];
% ax.XLim = [4, 20];
% ax.Title.Interpreter = 'LaTeX';
% ax.XLabel.Interpreter = 'LaTeX';
% ax.YLabel.Interpreter = 'LaTeX';
% ax.ZLabel.Interpreter = 'LaTeX';
% ax.TickLabelInterpreter = 'LaTeX';
% xlabel('$n_\perp$')
% ylabel('$\gamma_d/\gamma_b$')

%%

deltaperpvec=0.1:0.1:1.0;
nperpvec=[4,8,12,16,20];
Nnperp=length(nperpvec);

myBlue=[0.6,0.8,1];

figure
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

set(gca,'FontSize',45)
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
xlabel('$\delta_\perp/\lambda_0$')
ylabel('$n_\perp^4 \gamma_d/\gamma_b$')

%%

myBlue=[0.6,0.8,1];
figure
iterdperp=5;
omega0optvec=omega0opt_cell{iterdperp};
xiopt=xiopt_cell{iterdperp};

Nnperp=length(nperpvec);

mylegend=cell(Nnperp,1);
plot([0,1],[0,1],'w');
hold on
for iternperp=1:Nnperp
    L=(0.1*iterdperp*nperpvec(iternperp))^2*xiopt(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    omega0opt=omega0optvec(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    h=plot(L,omega0opt,'.-');
    h.MarkerSize=25;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    hold on
end

h=plot(L,sqrt(L/(2*pi)),'--r');
h.LineWidth=2;

set(gca,'FontSize',22)
ax = gca;
ax.XScale='log';
ax.YLim=[0,4];
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
grid on
xlabel('$L/\lambda_0$')
ylabel('$w_0$')

ax.Position=[0.1447,3*0.1501,0.7603,0.7749/2];

figure
for iternperp=1:Nnperp
    L=(0.1*iterdperp*nperpvec(iternperp))^2*xiopt(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    omega0opt=omega0optvec(iternperp,0.5*xiopt(iternperp,:).*(0.1*iterdperp*nperpvec(iternperp))^2>1);
    zR=pi*omega0opt.^2;
    h=plot(L,omega0opt.*sqrt(1+(0.5*L./zR).^2),'.-');
    h.MarkerSize=25;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    hold on
    
    h=plot(L,ones(size(L))*0.1*iterdperp*nperpvec(iternperp)/4.5,'--');
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
end
set(gca,'FontSize',22)
ax = gca;
ax.XScale='log';
ax.FontName = 'LaTeX';
ax.YLim=[0,4];
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
grid on
xlabel('$L/\lambda_0$')
ylabel('$w$')


h=plot(L,sqrt(L/(pi)),'--r');
h.LineWidth=2;
