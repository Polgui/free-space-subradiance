clearvars

load('Compress.mat')

nperp=20;

lambda0=1;

Nxi=20;

dperpvec=linspace(0.1,1.,10);
Ndperp=length(dperpvec);
        
gamma=1e1;

overlaps_vec=cell(Ndperp);
Delta_vec=cell(Ndperp);

for iterdperp=1:Ndperp
    
    dperp=dperpvec(iterdperp)
    
    xiopt=xiopt_cell{iterdperp};
    omega0opt=omega0opt_cell{iterdperp};

    Delta_temp=zeros(1,Nxi);
        
    for iterL=1:Nxi
        
        if 0.5*xiopt(5,iterL)*(nperp*dperp)^2 > 1

            L=0.5*xiopt(5,iterL)*(nperp*dperp)^2/lambda0+0.1;

            omega0=omega0opt(5,iterL);
            zR=pi*omega0^2/lambda0;

            xvec=setXpos(nperp,dperp);
            yvec=setYpos(nperp,dperp);
            zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

            gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
            [gammad,gammab,Deltalight,overlaps]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);
            
            Delta_temp(iterL)=Deltalight/gamma;
            
        else 
            
            Delta_temp(iterL)=NaN;
            
        end


    end
    Delta_vec{iterdperp}=Delta_temp;
end

%%

myBlue=[0.6,0.8,1];

figure
for iterdperp=3:Ndperp-1
    Delta_temp=Delta_vec{iterdperp};
    h=plot(xiopt(5,:),Delta_temp,'.-');    
    h.MarkerSize=25;
    h.LineWidth=2;
    h.Color=myBlue*(Ndperp-iterdperp+1)/Ndperp;
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

%%

