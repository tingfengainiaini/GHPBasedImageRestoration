function z   =  RDWT(x, wav, L)    % Redundant wavelet transform

% [t1,t2,lev]  =  mrdwt(x,wav,L);
% 
% [n  m]   =  size(t1);
% t1 = t1*2^(-L);
% for ll = 1:L
%     t2(:,(ll-1)*m*3+1:ll*m*3) = ...
%     t2(:,(ll-1)*m*3+1:ll*m*3)*2^(-ll);
% end
% 
% z            =  [t1 t2];

[t1,t2,lev]  =  mrdwt(x,wav,L);

[n  m]   =  size(t1);
% for ll = 1:L
%     t2(:,(ll-1)*m*3+1:ll*m*3) = ...
%     t2(:,(ll-1)*m*3+1:ll*m*3)*2^(-ll);
% end

z            =  [t1 t2];