function [Wu,Wy,b]=update_ddnet(y,Ud,Yd,Wu,Wy,b,sig,sigd,gam,nd)
ny=length(Yd)/nd;nu=length(Ud)/nd;
h=Wu*Ud+Wy*Yd+b;
ya=ddnet(Ud,Yd,Wu,Wy,b,sig);
Wub=Wu;
for k=1:ny
    it=1;
    for i=1:nd
        for j=1:nu           
            Wu(k,it)=Wub(k,it)-gam*(ya-y)*sigd(h(k))*Ud(it);
            it=it+1;
        end
    end
end
Wyb=Wy;
for k=1:ny
    it=1;
    for i=1:nd
        for j=1:ny           
            Wy(k,it)=Wyb(k,it)-gam*(ya-y)*sigd(h(k))*Yd(it);
            it=it+1;
        end
    end
end
bb=b;
for k=1:ny   
    b(k)=bb(k)-gam*(ya-y)*sigd(h(k));
end
end