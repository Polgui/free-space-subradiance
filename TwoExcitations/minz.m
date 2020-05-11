function val = minz(z,r2,n0,lambda0,zR)

val=abs(z-lambda0*atan(z/zR)/(2*pi)-n0*lambda0/4. + 0.5*r2./(z.*(1+(zR./z).^2)) );

end

