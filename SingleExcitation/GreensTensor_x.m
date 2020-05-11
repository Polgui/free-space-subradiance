function G = GreensTensor_x( lambda0,x,y,z )

r=sqrt(x.^2+y.^2);
R=sqrt(r.^2+z.^2);
k0R=2.*pi*R/lambda0;
G= exp(1j.*k0R).*(1.+1j./k0R-1./k0R.^2-(1.+3.*1j./k0R-3./k0R.^2).*x.*(x+1j*y)./(sqrt(2)*(r.^2+z.^2)))./(4.*pi.*R);

end