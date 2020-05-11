function xvec = setXpos( nperp,dperp )

xvec=zeros(1,nperp^2*2);

for iter=1:length(xvec)
    xvec(iter)=-0.5*dperp*(nperp-1)+mod(iter-1,nperp)*dperp;
end

end

