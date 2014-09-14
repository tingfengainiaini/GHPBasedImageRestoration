function z   = IRDWT(x, wav, L)    %Inverse redundant wavelet transform

% [rows  cols]   =  size(x);
% cols           =  (cols)/(L*3+1);
% 
% m   =  cols;
% 
% t1 = x(:,1:m)*2^L;
% for ll = 1:L
%      t2(:,(ll-1)*m*3+1:ll*m*3) = ...
%      x(:,m+(ll-1)*m*3+1:m+ll*m*3)*2^(ll);
% end
% 
% % t1             =  x(:, 1:cols);
% % t2             =  x(:, cols+1:end);
% 
% z              =  mirdwt(t1,t2,wav,L);


[rows  cols]   =  size(x);
cols           =  (cols)/(L*3+1);

m   =  cols;

t1 = x(:,1:m);
for ll = 1:L
     t2(:,(ll-1)*m*3+1:ll*m*3) = ...
     x(:,m+(ll-1)*m*3+1:m+ll*m*3);
end

z              =  mirdwt(t1,t2,wav,L);