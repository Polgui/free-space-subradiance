clearvars

load('Compress.mat')

dperp=.5;

xiopt=xiopt_cell{5};
omega0opt=omega0opt_cell{5};

lambda0=1;

Nxi=size(xiopt,2);

nperpvec=[4,8,12,16,20];
Nnperp=length(nperpvec);
        
gamma=1e1;

overlaps_vec=cell(Nnperp);
Delta_vec=cell(Nnperp);

for iternperp=1:Nnperp
    
    nperp=nperpvec(iternperp)
    
    overlap_temp=zeros(2*nperp^2,Nxi);
    Delta_temp=zeros(1,Nxi);
        
    for iterL=1:Nxi
        
        if 0.5*xiopt(iternperp,iterL)*(nperp*dperp)^2 > 1

            L=0.5*xiopt(iternperp,iterL)*(nperp*dperp)^2/lambda0+0.1;

            omega0=omega0opt(iternperp,iterL);
            zR=pi*omega0^2/lambda0;

            xvec=setXpos(nperp,dperp);
            yvec=setYpos(nperp,dperp);
            zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

            gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
            [gammad,gammab,Deltalight,overlaps]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);
            
            overlap_temp(:,iterL)=overlaps;
            Delta_temp(iterL)=Deltalight/gamma;
            
        else 
            
            overlap_temp(:,iterL)=ones(1,2*nperp^2)*NaN;
            Delta_temp(iterL)=NaN;
            
        end


    end
    overlaps_vec{iternperp}=overlap_temp;
    Delta_vec{iternperp}=Delta_temp;
end

%%

myBlue=[0.6,0.8,1];

figure
for iternperp=1:Nnperp
    overlap_temp=overlaps_vec{iternperp};
    for iterxi=1:Nxi
        overlap_sort=sort(abs(overlap_temp(:,iterxi)).^2);
        h=plot(xiopt(iternperp,iterxi),2*sum(overlap_sort(1:end-2)),'.-');    
        h.MarkerSize=45;
        h.LineWidth=2;
        h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
        hold on
        h=plot(xiopt(iternperp,iterxi),2*sum(overlap_sort(end-1:end)),'.-');    
        h.MarkerSize=45;
        h.LineWidth=2;
        h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    end
end
set(gca,'FontSize',45)
ax = gca;
ax.XLim = [1e-2, 1e1];
ax.YLim = [1e-4, 1e0];
ax.XTick=[1e-2,1e-1,1e0,1e1];
%ax.YTick=[1e-5,1e-4,1e-3,1e-2,1e-1,1e0];
ax.XScale='log';
ax.YScale='log';

ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
grid on
xlabel('$L\lambda_0/L_\perp^2$')
ylabel('Overlaps')

%%

myBlue=[0.6,0.8,1];

figure
for iternperp=1:Nnperp
    Delta_temp=Delta_vec{iternperp};
    h=plot(xiopt(iternperp,:),Delta_temp,'.-');    
    h.MarkerSize=25;
    h.LineWidth=2;
    h.Color=myBlue*(Nnperp-iternperp+1)/Nnperp;
    hold on
end
set(gca,'FontSize',22)
ax = gca;
ax.XLim = [1e-2, 1e1];
ax.XTick=[1e-2,1e-1,1e0,1e1];
ax.XScale='log';

ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
grid on
xlabel('$L\lambda_0/L_\perp^2$')
ylabel('$\Delta/\gamma$')