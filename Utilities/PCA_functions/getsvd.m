function [P, mx]=getsvd(X, nsig)

%X: MxN matrix (M dimensions, N trials)
%Y: Y=P*X
%P: the transform matrix
%V: the variance vector

% [M,N]=size(X);
% 
% mx   =  mean(X,2);
% mx2  =  repmat(mx,1,N);
% 
% X=X-mx2;
% 
% CovX=X*X'/(N-1);
% 
% CovX  =  CovX - diag(nsig^2*ones(size(CovX,1),1));
% ind = find(CovX<0); 
% CovX(ind) =0.0001; 
% 
% [P,V]=eig(CovX);
% 
% V=diag(V);
% [t,ind]=sort(-V);
% % V=V(ind);
% P=P(:,ind);
% P=P';
% % Y=P*X;
% 
% return;


mx   =  zeros(0);

[M,N]=size(X);

[P, S, V]   =  svd(X);
S=diag(S);
[t,ind]=sort(-S);
P=P(:,ind);
P=P';


% mx   =  mean(X,2);
% mx2  =  repmat(mx,1,N);
% 
% X=X-mx2;
% 
% CovX=X*X'/(N-1);
% 
% CovX  =  CovX - diag(nsig^2*ones(size(CovX,1),1));
% ind = find(CovX<0); 
% CovX(ind) =0.0001; 
% 
% [P,V]=eig(CovX);
% 
% V=diag(V);
% [t,ind]=sort(-V);
% % V=V(ind);
% P=P(:,ind);
% P=P';
% % Y=P*X;

return;
