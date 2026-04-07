function f=mpcerr(Uf,yd,Ud,Yd,Wu,Wy,b,sig)
ny=length(yd);nd=length(Yd)/ny;
nu=length(Ud)/nd;nf=length(Uf)/nu;
f=0;Yf=zeros(nf*ny,1);
for i=1:nf
    Ub=[Uf(1:(i-1)*nu);Ud(1:(nd-i+1)*nu)];
    Yb=[Yf(1:(i-1)*ny);Yd(1:(nd-i+1)*ny)];
    yf=ddnet(Ub,Yb,Wu,Wy,b,sig);
    Yf((i-1)*ny+1:i*ny)=yf;
    f=f+norm(yf-yd);
end
end