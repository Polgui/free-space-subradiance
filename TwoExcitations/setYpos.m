function yvec = setYpos( nperp,dperp )

yvec=zeros(1,nperp^2*2);

for iter=1:length(yvec)
    yvec(iter)=-0.5*dperp*(nperp-1)+mod(floor((iter-1)/nperp),nperp)*dperp;
end

end

