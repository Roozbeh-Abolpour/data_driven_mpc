function x=rmap(xh,kapa)
a=xh(end);
x=1/a*sqrt((kapa-1)/kapa)*xh(1:end-1);
end