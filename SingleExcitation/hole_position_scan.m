%% gammad/gammab
clearvars

lambda0=1;

dperp=0.5;
gamma=1e1;

nperp=12;
na=nperp^2;

Lvec=[5,7.5,10];

gammab_vec=zeros(length(Lvec),2*na);
gammad_vec=zeros(length(Lvec),2*na);
gvec_vec=zeros(length(Lvec),2*na);

gammad_free=zeros(length(Lvec),1);
gammab_free=zeros(length(Lvec),1);

for iterL=1:length(Lvec)
    
    L=Lvec(iterL)

    omega0=find_omega0opt(lambda0,nperp,dperp,L);
    zR=pi*omega0^2/lambda0;

    xvec=setXpos(nperp,dperp);
    yvec=setYpos(nperp,dperp);
    zvec=setZpos(nperp,L,lambda0,zR,xvec,yvec);

    gvec=Hopping(xvec,yvec,zvec,lambda0,omega0);
    [gammad_free(iterL),gammab_free(iterL)]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec);

    gvec_vec(iterL,:)=gvec;

    for iterhole=1:2*na

        Hnh=zeros(length(xvec),length(xvec));

        for j=1:length(xvec)
            Hnh(j,j)=Hnh(j,j)-(1-(j==iterhole))*1j*gamma/2.;
            for k=j+1:length(xvec)
               Jhop=(1-(j==iterhole))*(1-(k==iterhole))*1.5*gamma*lambda0*GreensTensor(lambda0,sqrt((xvec(k)-xvec(j))^2+(yvec(k)-yvec(j))^2),zvec(k)-zvec(j));
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

        Deltalight=real(maxoverlaps(indexgammad));

        Deltabright=real(maxoverlaps(3-indexgammad));

        overlaps=overlaps(index);

        gammad_vec(iterL,iterhole)=gammad;
        gammab_vec(iterL,iterhole)=gammab;
    end

end

%% Plots

myBlue=[0.6,0.8,1];
myRed=[.8,0.1,0.2];

figure
semilogy([0 .1],[0 .1],'w');
hold on
for iterL=length(Lvec):-1:1
    [gvecsort,index]=sort(abs(gvec_vec(iterL,:)).^2);
    ratiosort=gammad_vec(iterL,index)./gammab_vec(iterL,index);
    
    h=semilogy(ratiosort,'.');
    h.Color=myBlue*(iterL)/(length(Lvec));
    h.MarkerSize=20;
    hold on
end
mylegend{1}='$L/\lambda_0$';
mylegend{2}='$10$';
mylegend{3}='$15$';
mylegend{4}='$20$';
hl=legend(mylegend);
hl.Interpreter='LaTeX';
hl.Position=[0.753 0.1763 0.1172 0.3549];
for iterL=1:length(Lvec)
    [gvecsort,index]=sort(abs(gvec_vec(iterL,:)).^2);
    ratiosort=gammad_vec(iterL,index)./gammab_vec(iterL,index);

    h=semilogy((gamma*gvecsort+gammad_free(iterL)*(1-2*gvecsort))./(gamma*gvecsort+gammab_free(iterL)*(1-2*gvecsort)),'.');
    h.MarkerSize=15;
    h.Color=myRed;
    hold on
end

sortindex=7:-1:1;
chH = get(gca,'Children');
set(gca,'Children',chH(sortindex))

set(gca,'FontSize',22)
ax = gca;
ax.FontName = 'LaTeX';
grid on
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.ZLabel.Interpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
xlabel('$j$')
ylabel('$\gamma_d/\gamma_b$')