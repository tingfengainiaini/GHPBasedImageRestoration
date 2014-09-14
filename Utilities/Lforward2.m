function X=Lforward2(P1,P2)

[m2,n2]=size(P1);
[m1,n1]=size(P2);

if (n2~=n1+1)
    error('dimensions are not consistent')
end
if(m1~=m2+1)
    error('dimensions are not consistent')
end

m=m2+1;
n=n2;

X=zeros(m,n);
X(1:m-1,:)=P1;
X(:,1:n-1)=X(:,1:n-1)+P2;
X(2:m,:)=X(2:m,:)-P1;
X(:,2:n)=X(:,2:n)-P2;


