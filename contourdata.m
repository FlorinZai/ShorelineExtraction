function s=contourdata(c)
tol=1e-12;
k=1;
col=1;
while col<size(c,2);
    s(k).level=c(1,col);
    s(k).num=c(2,col);
    idx=col+1:col+c(2,col);
    s(k).xdata=c(1,idx);
    s(k).ydata=c(2,idx);
    s(k).isopen=abs(diff(c(1,idx([1 end]))))>tol||abs(diff(c(2,idx([1 end]))))>tol;
    k=k+1;
    col=col+c(2,col)+1;
end
    


end

