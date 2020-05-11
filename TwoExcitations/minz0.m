function val = minz0(z,lambda0,zR)

val=abs(mod((z-lambda0*atan(z/zR)/(2*pi)),lambda0/4.));

end

