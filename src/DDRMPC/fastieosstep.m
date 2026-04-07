function xh=fastieosstep(Qhs,Ah,xh)
conv_tol=1e-4;eq_tol=1e-4;
while 1==1
    d=0;
    for i=1:length(Qhs)
        d=max(d,xh'*Qhs{i}*xh);
    end
    xh=xh/sqrt(d);
    G=[];indcs=[];
    for i=1:length(Qhs)
        if xh'*Qhs{i}*xh>=1-eq_tol
            G=[G Qhs{i}*xh];
            indcs=[indcs i];
        end
    end
    l=fastquad(G'*G+1e-6*eye(size(G,2)));
    xt=G*l;
    if xt'*xt>1/(xh'*xh)
        a=0.5*(xt'*xt+1/(xh'*xh));
        p=a*xh-xt;
    else
        a=0.5*(xt'*xt+1/(xh'*xh));
        p=-a*xh+xt;
    end
    s=10;p=p/norm(p);
    for i=1:length(Qhs)
        rs=roots([p'*Qhs{i}*p 2*xh'*Qhs{i}*p xh'*Qhs{i}*xh-1]);
        s=min(s,max(rs));
    end

    if xh'*xh>=1-conv_tol||s<=conv_tol
        return
    end  
    xh=xh+s*p;

end
end