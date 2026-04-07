function vh=ieosstep(Qhs,Ah,vhb)
m=length(Qhs);n=length(Qhs{1});
itmax=100;it=1;
while it<=itmax
    yalmip('clear')
    vh=sdpvar(n,1);
    if ~isempty(Ah)
        C=[Ah*vh<=0];
    else
        C=[];
    end
    for i=1:m
        C=[C,vh'*Qhs{i}*vh<=1];
    end
    C=[C,vh(end-1)<=vh(end)*vhb(end-1)/vhb(end),vh(end)>=1e-4];
    op=sdpsettings;op.solver='SDPT3';op.verbose=0;
    res=optimize(C,-vhb'*vh,op);
    vh=value(vh);    
    
    if vh'*vh>=1-1e-5||norm(vh-vhb)<=1e-4
        break
    end
    vhb=vh;
    it=it+1;
end
end