function u=nnddmpc(u,y,yref,bounds,orders,k)
global Ud Yd Wu Wy Ufb b;
yd=yref;gam=1e-2;
umin=bounds.umin;umax=bounds.umax;
ymin=bounds.ymin;ymax=bounds.ymax;
udmax=bounds.dumax;
nd=orders.nd;nf=orders.nf;
nu=1;ny=1;
sig=@(x) 20*atan(0.1*x)/pi;
x=sym('x');
str=sprintf('sigd=@(x) %s;',diff(sig(x)));
eval(str)
if k==0
    Ufb=[];
    Wu=1-2*rand(1,nu*nd);
    Wy=1-2*rand(1,ny*nd);
    b=1-2*rand;
    Ud=zeros(nd*nu,1);
    Yd=zeros(nd*ny,1);
end
Ud=[u;Ud(1:end-nu)];
Yd=[y;Yd(1:end-ny)];
if k>nd
    fptr=@(Uf) mpcerr(Uf,yd,Ud,Yd,Wu,Wy,b,sig);
    cptr=@(Uf) mpccon(Uf,yd,Ud,Yd,Wu,Wy,b,sig,umin,umax,ymin,ymax,udmax);
    if isempty(Ufb)
        Uf0=zeros(nf*nu,1);
    else
        Uf0=Ufb;
    end
    try
    Uf=fmincon(fptr,Uf0,[],[],[],[],[],[],cptr,optimset('disp','none'));
    u=Uf(1:nu);
    Ufb=Uf;
    catch
        u=0.1*(umin+rand*(umax-umin));
    end
    for i=1:nu
        if u(i)<umin(i)
            u(i)=umin(i);
        end
        if u(i)>umax(i)
            u(i)=umax(i);
        end
    end
else
    u=0.1*(umin+rand*(umax-umin));
end
[Wu,Wy,b]=update_ddnet(y,Ud,Yd,Wu,Wy,b,sig,sigd,gam,nd);
end