function X1=skel_HMT(X, B)

if size(size(X),2)==2
    H=ones(3,3);
end

if size(size(X),2)==3
    H=ones(3,3,3);
end

SE1=B;
SE2=H-SE1;
X1=imerode(X,SE1) & imerode(~X,SE2);
