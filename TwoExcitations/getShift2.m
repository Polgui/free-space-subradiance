function [gammad2,eigvec2,sortoverlaps,Deltalight2]=getShift2(gvec,lambda0,gamma,xvec,yvec,zvec)
na=length(xvec)/2;

dim=na*(na*2-1);

Hnh=zeros(dim,dim);

for j=1:2*na
   for k=j+1:2*na
        index= 2*na*(j-1)-j*(j-1)/2+(k-j-1);
        Hnh(index+1,index+1)=Hnh(index+1,index+1)-1j*gamma;

        for l=k+1:2*na
            index1= 2*na*(j-1)-j*(j-1)/2+(l-j-1);
            index2= 2*na*(k-1)-k*(k-1)/2+(l-k-1);

            Jjk=1.5*gamma*lambda0*GreensTensor(lambda0,sqrt((xvec(k)-xvec(j))^2+(yvec(k)-yvec(j))^2),zvec(k)-zvec(j));
            Jkl=1.5*gamma*lambda0*GreensTensor(lambda0,sqrt((xvec(l)-xvec(k))^2+(yvec(l)-yvec(k))^2),zvec(l)-zvec(k));
            Jjl=1.5*gamma*lambda0*GreensTensor(lambda0,sqrt((xvec(l)-xvec(j))^2+(yvec(l)-yvec(j))^2),zvec(l)-zvec(j));

            Hnh(index+1,index1+1)=Hnh(index+1,index1+1)-Jkl;
            Hnh(index1+1,index+1)=Hnh(index1+1,index+1)-Jkl;

            
            Hnh(index+1,index2+1)=Hnh(index+1,index2+1)-Jjl;
            Hnh(index2+1,index+1)=Hnh(index2+1,index+1)-Jjl;


            Hnh(index1+1,index2+1)=Hnh(index1+1,index2+1)-Jjk;
            Hnh(index2+1,index1+1)=Hnh(index2+1,index1+1)-Jjk;

        end
   end
end

[u,v]=eig(Hnh);
lam=diag(v);

gvecdark=gvec.*(-1).^(zvec>=0);

gdark=zeros(dim,1);
for j=1:2*na
   for k=j+1:2*na
        index= 2*na*(j-1)-j*(j-1)/2+(k-j-1);
        gdark(index+1)=gvecdark(j)*gvecdark(k)+gvecdark(j)*gvecdark(k);
    end
end
gdark=gdark/norm(gdark);

overlaps=gdark'*u;
[sortoverlaps,index]=sort(abs(overlaps).^2);

gammad2=-2*imag(lam(index(end)));

Deltalight2=real(lam(index(end)));

eigvec2=u(:,index(end));

end

