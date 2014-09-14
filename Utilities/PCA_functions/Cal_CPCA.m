function [CPCA_D mX] =  Cal_CPCA( im, par )

b        =  par.win;
t        =  floor( b/2 );
T        =  par.T;
st       =  par.step;
[h w ch] =  size(im);
num      =  ceil(h/st)*ceil(w/st);

CPCA_D   =  zeros(b*b*b*b*ch*ch, num, 'single');
mX       =  zeros(b*b*ch, num, 'single');
cnt      =  1;

im       =  single(im);
for  row  =  1:st:h

    if  (row+b-1) > h
        row   =  h-b+1;
    end       
    
    for  col  =  1:st:w
                
        if  (col+b-1) > w
            col   =  w-b+1;
        end                      

        cb      =  im(row:row+b-1, col:col+b-1, :);
        rmin    =  max( row+t-T, 1);
        rmax    =  min( row+t+T, h);
        cmin    =  max( col+t-T, 1);
        cmax    =  min( col+t+T, w);
        Blk     =  im(rmin:rmax, cmin:cmax, :);
        
        [n, m]  =  size(Blk(:,:,1));
        N       =  n-b+1;
        M       =  m-b+1;
        L       =  N*M;
        X       =  zeros(b*b*ch,L);

        cb2(:,:,1)  =  cb(:,:,1)';
        cb2(:,:,2)  =  cb(:,:,2)';
        cb2(:,:,3)  =  cb(:,:,3)';
        cb          =  cb2;
        cb          =  cb(:);
        cb          =  repmat(cb,1,L);

        w2  =  b*b;
        k   =  0;
        X   =  zeros(b*b*ch,L);
        for i  = 1:b
            for j  = 1:b
                k     =  k+1;
                blk   =  Blk(i:n-b+i,j:m-b+j,:);
                
                blk1  =  blk(:,:,1);
                blk2  =  blk(:,:,2);
                blk3  =  blk(:,:,3);
                
                blk1         =  blk1(:);
                X(k,:)       =  blk1';
                blk2         =  blk2(:);
                X(k+w2,:)    =  blk2';
                blk3         =  blk3(:);
                X(k+w2*2,:)  =  blk3';
                
%                 blk  =  blk(:);
%                 X(k,:) =  blk';
            end
        end

        E    =  abs(X-cb).^2;
        mE   =  mean(E);
        [val, ind]   =  sort(mE);
        
        num  =  min( par.num, L );
        X    =  X(:,ind(1:num));
        
        [P, mx]   =  getpca(X);
        
        CPCA_D(:,cnt)   =  single(P(:));
        mX(:,cnt)       =  single(mx);
        cnt             =  cnt + 1;
    end
end
        