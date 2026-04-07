function [F,Fe]=mpccon(Uf,yd,Ud,Yd,Wu,Wy,b,sig,umin,umax,ymin,ymax,udmax)
ny=length(yd);nd=length(Yd)/ny;
nu=length(Ud)/nd;nf=length(Uf)/nu;
F=[];Fe=[];Yf=zeros(nf*ny,1);
for i=1:nf    
    Ub=[Uf(1:(i-1)*nu);Ud(1:(nd-i+1)*nu)];
    Yb=[Yf(1:(i-1)*ny);Yd(1:(nd-i+1)*ny)];
    yf=ddnet(Ub,Yb,Wu,Wy,b,sig);
    Yf((i-1)*ny+1:i*ny)=yf;
    for j=1:ny
        F=[F;Yf((i-1)*ny+1:i*ny)-ymax(j)];
        F=[F;-Yf((i-1)*ny+1:i*ny)+ymin(j)];
    end
    for j=1:nu
        F=[F;Uf((i-1)*nu+1:i*nu)-umax(j)];
        F=[F;-Uf((i-1)*nu+1:i*nu)+umin(j)];
    end
end
for j=1:nu
    F=[F;Uf(j)-Ub(j)-udmax(j)];
    F=[F;-Uf(j)+Ub(j)-udmax(j)];
end
for i=2:nf-1
    for j=1:nu
        F=[F;Uf((i-1)*nu+j)-Uf(i*nu+j)-udmax(j)];
        F=[F;Uf(i*nu+j)-Uf((i-1)*nu+j)-udmax(j)];
    end
end
end