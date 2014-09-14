function [PCA_D arr_idx mX] = Clustering_LPCA( im, opts )

[h w ch]  = size( im );

if ch==3
   y_im      =   rgb2ycbcr( uint8(im) );
   y_im      =   double( y_im(:,:,1) );
end

cls      =  opts.cls;
b        =  opts.win;
hpf      =  fspecial('gauss', 7, 1);
hp_im    =  Blur( 'fwd', y_im, hpf );
hp_im    =  y_im - hp_im;


% Clustering the lum component of the input image
N       =  h-b+1;
M       =  w-b+1;
L       =  N*M;
X       =  zeros(b*b,L, 'single');
k       =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        blk  =  hp_im(i:h-b+i,j:w-b+j);
        blk  =  blk(:);
        X(k,:) =  blk';
    end
end

P      =  randperm(L);
P2     =  P(1:cls);
vec    =  X(:,P2(1:end));
L2     =  ceil(L*opts.rat);
if L2<(cls*opts.num)
    L2   = L;
end
P      =  P(1:L2);

% Do the clustering
for i = 1:opts.itn
    mse       =  0;
    cent      =  zeros(b*b, cls);
    cnt       =  zeros(1,cls, 'single');
    arr_idx   =  zeros(round(L/20), cls, 'single');
    for j = 1:L2
        v     =  X(:,P(j));
        cb    =  repmat(v,1,cls);
        dis   =  sum((vec-cb).^2);
        [val ind]      =  min( dis );
        cent(:, ind)   =  cent(:,ind) + v;
        cnt( ind )     =  cnt( ind ) + 1;
        arr_idx(cnt(ind), ind)  =  P(j);
        
        mse  =  mse + val;
    end
    
    if i==ceil(opts.itn-1)
        [val ind]  =  min( cnt );       % Remove these classes with little samples
        while val<opts.num
            vec(:,ind)  = [];
            cent(:,ind) = [];
            cls      =  cls - 1;
            arr      =  arr_idx(:, ind);
            arr_idx(:,ind)  =  [];
            cnt(ind)  =  [];
            for p = 1:val
                j    =  arr(p);
                v    =  X(:, j);
        
                cb      =  repmat(v,1,cls);
                dis     =  sum((vec-cb).^2);
                [val k]     =  min( dis );
                cent(:,k)   =  cent(:,k) + v;
                cnt( k )    =  cnt( k ) + 1;
                arr_idx(cnt(k),k)  =  j;
            end
    
            [val  ind]   =  min(cnt);
        end        
    end
    
    cnt   =  cnt + eps;
    wei   =  repmat(cnt, b*b, 1);
    vec   =  cent./wei;    
    mse   =  mse/L2/b/b;
%     disp( sprintf('clustering %d th loop, mse = %f', i, mse) );
    
end


arr_idx   =  zeros(round(L/20), cls, 'single');
cnt       =  zeros(1,cls, 'single');
for i = 1:L
    
    v    =  X(:,i);
    cb   =  repmat(v,1,cls);
    dis  =  sum( (vec-cb).^2 );
    [val ind]     =  min( dis );
    
    cnt( ind )    =  cnt( ind ) + 1;
    arr_idx(cnt(ind), ind)  =  i;
end


% Remove these classes with little samples
[val ind]  =  min( cnt );
while val<opts.num
    vec(:,ind)  = [];
    cls      =  cls - 1;
    arr      =  arr_idx(:, ind);
    arr_idx(:,ind)  =  [];
    cnt(ind)  =  [];
    for i = 1:val
        j    =  arr(i);
        v    =  X(:, j);
        
        cb      =  repmat(v,1,cls);
        dis     =  sum((vec-cb).^2);
        [val k]     =  min( dis );
        cnt( k )    =  cnt( k ) + 1;
        arr_idx(cnt(k),k)  =  j;
    end
    
    [val  ind]   =  min(cnt);
end


% Train the color/gray PCAs
k   =  0;
X   =  zeros(b*b*ch,L, 'single');
if ch==1
    for i  = 1:b
        for j  = 1:b
            k    =  k+1;
            blk  =  im(i:h-b+i,j:w-b+j);
            blk  =  blk(:);
            X(k,:) =  blk';
        end
    end
else
    % RGB->YUV 
    y_im   =  double( rgb2ycbcr( uint8(im) ) );
    
    for i  = 1:b
        for j  = 1:b
            k    =  k+1;
            blk  =  y_im(i:h-b+i,j:w-b+j, :);
            
            blk1  =  blk(:,:,1);
            blk2  =  blk(:,:,2);
            blk3  =  blk(:,:,3);
                
            blk1         =  blk1(:);
            X(k,:)       =  blk1';
            blk2         =  blk2(:);
            X(k+b*b,:)   =  blk2';
            blk3         =  blk3(:);
            X(k+b*b*2,:) =  blk3';
        end
    end
    
end

PCA_D   =  zeros(b*b*b*b*ch*ch, cls, 'single');
mX      =  zeros(b*b*ch, cls, 'single');
for  i = 1:cls
    ind   =  arr_idx(:, i)';
    ind   =  ind( 1:cnt(i) );
    X1    =  X(:, ind(1:end));

    [P, mx]   =  getpca(X1);
        
    PCA_D(:,i)    =  P(:);
    mX(:,i)       =  mx;    
end
    