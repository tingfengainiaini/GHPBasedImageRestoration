function  imout  =  Soft_LPCA( im, PCA_D, arr_idx, mX, opts )

[h w ch]   =  size(im);
b     =  opts.win;
tau   =  opts.tau;
cls   =  size(PCA_D, 2);
k     =  0;

N     =  h-b+1;
M     =  w-b+1;
L     =  N*M;
X     =  zeros(b*b*ch,L,'single');
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

for i = 1:cls
    ind        =  arr_idx(:, i)';
    [val idx]  =  find(ind);
    ind   =  ind( idx );
    cnt   =  size(ind,2);
    X1    =  X(:, ind(1:end));
    m1    =  repmat(mX(:,i), 1, cnt);    

    P     =  reshape( PCA_D(:, i), (b^2)*ch, (b^2)*ch );
    Y1    =  soft(P*(X1-m1), tau);
    X1    =  P'*Y1 + m1;
    X(:, ind(1:end))  =  X1;    
end

% Output the processed image
im_out   =  zeros( size(im) );
im_wei   =  zeros(h,w);
k        =  0;
if ch==1
    for i  = 1:b
        for j  = 1:b
            k    =  k+1;
            blk  =  X(k,:)';
            blk  =  reshape(blk, [N M]);
            im_out(i:h-b+i,j:w-b+j)  =  im_out(i:h-b+i,j:w-b+j) + blk;
            im_wei(i:h-b+i,j:w-b+j)  =  im_wei(i:h-b+i,j:w-b+j) + 1;
        end
    end
else
    for i  = 1:b
        for j  = 1:b
            k     =  k+1;
            blk1  =  X(k,:)';
            blk2  =  X(k+b*b,:)';
            blk3  =  X(k+b*b*2,:)';
            blk1  =  reshape(blk1, [N M]);
            blk2  =  reshape(blk2, [N M]);
            blk3  =  reshape(blk3, [N M]);
            blk(:,:,1)  =  blk1;
            blk(:,:,2)  =  blk2;
            blk(:,:,3)  =  blk3;
            im_out(i:h-b+i,j:w-b+j,:)  =  im_out(i:h-b+i,j:w-b+j,:) + blk;
            im_wei(i:h-b+i,j:w-b+j)    =  im_wei(i:h-b+i,j:w-b+j) + 1;
        end
    end 
end

if ch==1
    imout  =  im_out./im_wei;
else
    im_out(:,:,1)   =  im_out(:,:,1)./im_wei;
    im_out(:,:,2)   =  im_out(:,:,2)./im_wei;
    im_out(:,:,3)   =  im_out(:,:,3)./im_wei;
    % YUV->RGB
    imout   =  double(ycbcr2rgb( uint8(im_out) ));
end

