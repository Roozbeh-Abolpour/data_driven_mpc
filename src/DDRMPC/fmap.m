function xh=fmap(x,kapa,Qs)
w=0;
for i=1:length(Qs)
    w=max(w,[x;1]'*Qs{i}*[x;1]-1);
end
a=sqrt(kapa-1)/sqrt(kapa*(x'*x)+kapa+w);
xh=[sqrt(kapa)/sqrt(kapa-1)*a*x;a];
end