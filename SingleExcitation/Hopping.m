function [Hopvec,GaussianBeam] = Hopping(xvec,yvec,zvec,lambda0,omega0 )

r2=xvec.^2+yvec.^2;
zR=pi*omega0^2/lambda0;
omegaz=omega0*sqrt(1+(zvec/zR).^2);
Rz=zvec.*(1+(zR./zvec).^2);
psiz=atan(zvec/zR);

GaussianBeam=sqrt(2/pi)*(1./omegaz).*exp(-r2./(omegaz.^2)).*exp(-1j*psiz).*exp(2*pi*1j*(zvec+r2./(2*Rz))/lambda0);

Hopvec=GaussianBeam/norm(GaussianBeam);
end

