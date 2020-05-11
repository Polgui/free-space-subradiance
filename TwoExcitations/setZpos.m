function [zvec,n0] = setZpos( nperp,L,lambda0,zR,xvec,yvec,dL )

if nargin < 7
    dL=0;
end

zvec=zeros(1,nperp^2*2);

zparse=(L-lambda0/4:lambda0/10000:L);
z0vec=minz0(zparse,lambda0,zR);
[~,argmin]=min(z0vec);
z0=zparse(argmin);

n0=(z0-lambda0*atan(z0/zR)/(2*pi))/(lambda0/4.);

for iter=1:2*nperp^2
    if iter<=nperp^2
        r2=xvec(iter)^2+yvec(iter)^2;
        
        zparse=(z0-3*lambda0:lambda0/20000:z0+lambda0);
        z0vec=minz(zparse,r2,n0,lambda0,zR);
        [~,argmin]=min(z0vec);

        zvec(iter)=-zparse(argmin)-dL;
    else
        zvec(iter)=-zvec(iter-nperp^2);
    end
            
end

end