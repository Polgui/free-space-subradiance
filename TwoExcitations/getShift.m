function [gammad,gammab,Deltalight,overlaps]=getShift(gvec,lambda0,gamma,xvec,yvec,zvec)

Hnh=zeros(length(xvec),length(xvec));
        
for j=1:length(xvec)
    Hnh(j,j)=Hnh(j,j)-1j*gamma/2.;
    for k=j+1:length(xvec)
       Jhop=1.5*gamma*lambda0*GreensTensor(lambda0,sqrt((xvec(k)-xvec(j))^2+(yvec(k)-yvec(j))^2),zvec(k)-zvec(j));
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

gammab=max(-2*imag(maxoverlaps));
[gammad,indexgammad]=min(-2*imag(maxoverlaps));

Deltalight=real(maxoverlaps(indexgammad));

end

