function x=fastquad(Q)
n=length(Q);e=ones(n,1);I=eye(n);
Qi=Q\eye(n);x=Qi*e/(e'*Qi*e);
for i=1:n         
    [o,k]=min(x);
    if o>=0        
        return
    end
    Ek=I;ek=I(:,k);
    Ek(:,k)=[];
    Qi=Ek*(Ek'*Qi*Ek-Ek'*Qi*(ek*ek')*Qi*Ek/Qi(k,k))*Ek';
    x=Qi*e/(e'*Qi*e);
end
end