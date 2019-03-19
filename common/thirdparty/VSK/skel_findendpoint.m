function [endpoint,X1,X2]=skel_findendpoint(A)

s=zeros(3,3,3);
s(2,2,2)=1;

rank=1;
for i=1:3
    for j=1:3
        for k=1:3
            B(:,:,:,rank)=s;
            if ~(i==j && i==k && j==k && i==2)
                B(i,j,k,rank)=1;
                rank=rank+1;
            end
        end
    end
end
     
SE=B(:,:,:,i);
X1=A-skel_HMT(A,SE);

X2=skel_HMT(X1,B(:,:,:,1));

for i=2:rank-1
    X2=X2 | skel_HMT(X1,B(:,:,:,i));
end

d1=reshape(X2,prod(size(X2)),1);
[indx3,indy3,indz3]=ind2sub(size(X2),find(d1==1));
endpoint=[indx3,indy3,indz3];

