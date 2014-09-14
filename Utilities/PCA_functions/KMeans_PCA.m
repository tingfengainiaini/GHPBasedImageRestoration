function  Dict   =  KMeans_PCA( im, par, cls_num )
b         =   par.win;
psf       =   fspecial('gauss', par.win+2, par.sigma);
[Y, X]    =   Get_patches(im, b, psf, 1, 1);
for i=2:2:6  % 1:1:6
   [Ys, Xs]   =   Get_patches(im, b, psf, 1, 0.8^i);
   Y          =   [Y Ys];
   X          =   [X Xs];
end

% Compute PCA for the smooth patches
delta       =   sqrt(par.nSig^2+16);
v           =   sqrt( mean( Y.^2 ) );
[a, i0]     =   find( v<delta );
if length(i0)<par.m_num
    D0      =   dctmtx(par.win^2);
else
    D0      =   getpca( X(:, i0), par.nSig );
end
set         =   1:size(Y, 2);
set(i0)     =   [];
Y           =   Y(:,set);
X           =   X(:,set);

% Clustering
b2        =   size(Y, 1);
L         =   size(Y, 2);
itn       =   14;
m_num     =   par.m_num;
rand('seed',0);


% % time of matlab implementation 
% time0       =   clock;
% opts = statset('MaxIter', itn);
% [cidx, ctrs] = kmeans(Y', cls_num,  'Options',opts);
% disp(sprintf('Total elapsed time of matlab k-means   = %f sec\n', (etime(clock,time0)) ));

% n0     =   ceil(0.5*b2);
% P0     =   getpca(Y, par.nSig);
% Y2     =   P0*Y;
% Y      =   Y2(1:n0, :);


% time of my own implementation 
% time0       =   clock;
[cls_idx, vec, cls_num]   =  My_kmeans(Y, cls_num, itn, m_num);
vec   =  vec';
% disp(sprintf('Total elapsed time of my own k-means   = %f sec\n', (etime(clock,time0)) ));


% Compute the PCA basis for eacl cluster
[s_idx, seg]   =  Proc_cls_idx( cls_idx );
PCA_D          =  zeros(b^4, cls_num);
for  i  =  1 : length(seg)-1
   
    idx    =   s_idx(seg(i)+1:seg(i+1));    
    cls    =   cls_idx(idx(1));    
    X1     =   X(:, idx);

    [P, mx]   =  getpca(X1, par.nSig);    %getpca(X1, par.nSig);    
    PCA_D(:,cls)    =  P(:);
end


% Assign each patch a PCA dictionary
[Y, X]      =   Get_patches(im, b, psf, par.step, 1);
cls_idx     =   zeros(size(X, 2), 1);

v           =   sqrt( mean( Y.^2 ) );
[a, ind]    =   find( v<delta );

% Y           =   P0*Y;
% Y           =   Y(1:n0,:);

% cy          =   D0*X(:, ind);
% cy          =   mean(cy.^2, 2);
% e_nsig      =   sqrt(mean( cy(9:end),1));

set         =   1:size(X, 2);
set(ind)    =   [];
L           =   size(set,2);
vec         =   vec';
b2          =   size(Y, 1);

for j = 1 : L
    dis   =   (vec(:, 1) -  Y(1, set(j))).^2;
    for i = 2 : b2
        dis  =  dis + (vec(:, i)-Y(i, set(j))).^2;
    end
    [val ind]      =   min( dis );
    cls_idx( set(j) )   =   ind;
end

[s_idx, seg]   =  Proc_cls_idx( cls_idx );

Dict.PCA_D       =   PCA_D;
Dict.cls_idx     =   cls_idx;
Dict.s_idx       =   s_idx;
Dict.seg         =   seg;
Dict.D0          =   D0;


%--------------------------------------------------
function  [Py  Px]  =  Get_patches( im, b, psf, s, scale )
im        =  imresize( im, scale, 'bilinear' );

[h w ch]  =  size(im);
ws        =  floor( size(psf,1)/2 );

if  ch==3
    lrim      =  rgb2ycbcr( uint8(im) );
    im        =  double( lrim(:,:,1));    
end

lp_im     =  conv2( psf, im );
lp_im     =  lp_im(ws+1:h+ws, ws+1:w+ws);
hp_im     =  im - lp_im;

N         =  h-b+1;
M         =  w-b+1;
% s         =  1;
r         =  [1:s:N];
r         =  [r r(end)+1:N];
c         =  [1:s:M];
c         =  [c c(end)+1:M];
L         =  length(r)*length(c);
Py        =  zeros(b*b, L, 'single');
Px        =  zeros(b*b, L, 'single');

k    =  0;
for i  = 1:b
    for j  = 1:b
        k       =  k+1;
        blk     =  hp_im(r-1+i,c-1+j);
        Py(k,:) =  blk(:)';
        
        blk     =  im(r-1+i,c-1+j);
        Px(k,:) =  blk(:)';        
    end
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
function   [cls_idx,vec,cls_num]  =  My_kmeans(Y, cls_num, itn, m_num)
Y         =   Y';
[L b2]    =   size(Y);
P         =   randperm(L);
P2        =   P(1:cls_num);
vec       =   Y(P2(1:end), :);

for i = 1 : itn
    
    mse       =  0;
    cnt       =  zeros(1, cls_num);    
    
    v_dis    =   zeros(L, cls_num);
    for  k = 1 : cls_num
        v_dis(:, k) = (Y(:,1) - vec(k,1)).^2;
        for c = 2:b2
            v_dis(:,k) =  v_dis(:,k) + (Y(:,c) - vec(k,c)).^2;
        end
    end

    [val cls_idx]     =   min(v_dis, [], 2);
    mse               =   mse + sum( val );
    
    [s_idx, seg]   =  Proc_cls_idx( cls_idx );
    for  k  =  1 : length(seg)-1
        idx    =   s_idx(seg(k)+1:seg(k+1));    
        cls    =   cls_idx(idx(1));    
        vec(cls,:)    =   mean(Y(idx, :));
        cnt(cls)      =   length(idx);
    end        
    
    if (i==itn-2)
        [val ind]  =  min( cnt );       % Remove these classes with little samples        
        while (val<m_num) && (cls_num>=40)
            vec(ind, :)    =  [];
            cls_num       =  cls_num - 1;
            cnt(ind)      =  [];

            [val  ind]    =  min(cnt);
        end        
    end
    
    mse   =  mse/L/b2;
%     disp( sprintf('clustering %d th loop, mse = %f, cls_num = %d', i, mse, cls_num) );    
end


% %------------------------------------------------------------------------
% function   [cls_idx,vec]  =  My_kmeans(Y, cls_num, itn)
% Y         =   Y';
% 
% [L b2]    =   size(Y);
% P         =   randperm(L);
% P2        =   P(1:cls_num);
% vec       =   Y(P2(1:end), :);
% T         =   60000;
% 
% for i = 1 : itn
%     
%     mse       =  0;
%     cnt       =  zeros(1, cls_num);
%     cls_idx   =  zeros(L, 1);
%     
%     for  j  =  1 : T : L        
%         js     =   min(j+T-1, L);
%         n      =   js-j+1;
%         
%         v_dis    =   zeros(n, cls_num, 'single');
%         for  k  =  1 : cls_num            
%             cv     =   vec(k,:);
%             cv     =   cv(ones(n,1), :);
% %             cv     =   repmat( vec(:,k), 1, n );    
%             v_dis(:, k)    =   sum( (Y(j:js, :)-cv).^2, 2 );            
%         end        
%         [val ind]         =   min(v_dis');
%         cls_idx(j:js,1)   =   ind;
%         mse               =   mse + sum( val );
%     end
%     
%     [s_idx, seg]   =  Proc_cls_idx( cls_idx );
%     for  k  =  1 : length(seg)-1
%         idx    =   s_idx(seg(k)+1:seg(k+1));    
%         cls    =   cls_idx(idx(1));    
%         vec(cls,:)    =   mean(Y(idx, :));
% %         cnt(cls)      =   length(idx);
%     end          
%     
% %     if (i==4 || i==itn-2)
% %         [val ind]  =  min( cnt );       % Remove these classes with little samples        
% %         while (val<m_num) && (cls_num>=40)
% %             vec(:,ind)    =  [];
% %             cls_num       =  cls_num - 1;
% %             cnt(ind)      =  [];
% % 
% %             [val  ind]    =  min(cnt);
% %         end        
% %     end
%     
%     mse   =  mse/L/b2;
%     disp( sprintf('clustering %d th loop, mse = %f', i, mse) );    
% end