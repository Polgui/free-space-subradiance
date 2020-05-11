function G = GreensTensor( lambda0,r,z )

R=sqrt(r.^2+z.^2);
k0R=2.*pi*R/lambda0;
G= exp(1j.*k0R).*(1.+1j./k0R-1./k0R.^2-(1.+3.*1j./k0R-3./k0R.^2).*r.^2./(2.*r.^2+2.*z.^2))./(4.*pi.*R);

end