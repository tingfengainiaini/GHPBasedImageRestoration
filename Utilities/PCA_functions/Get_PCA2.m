function  PCA_idx =  Get_PCA2( im, par, codewords )

[h  w]     =  size(im);
cls_num    =  size( codewords, 2 );
LP_filter  =  fspecial('gaussian', 7, par.sigma);
lp_im      =  conv2( LP_filter, im );
lp_im      =  lp_im(4:h+3, 4:w+3);
hp_im      =  im - lp_im;

b       =  par.win;
t       =  floor( b/2 );
st      =  par.step;
len     =  b*b;


PCA_idx  =  zeros(h*w,1);
cnt      =  1;

%-------------
[P, mx]     =  getpca(codewords);
mx          =  mean(codewords,2);
mx2         =  repmat(mx,1,cls_num);
codewords   =  codewords-mx2;
num         =  25;
P           =  P(1:num, :);
codewords   =  P(1:num, :)*codewords;


for  row  =  1:st:h

    if  (row+b-1) > h
        row   =  h-b+1;
    end       
    
    for  col  =  1:st:w
                
        if  (col+b-1) > w
            col   =  w-b+1;
        end                      

        cb             =   hp_im(row:row+b-1, col:col+b-1);
        mm             =   mean(mean(cb));
        vec            =   P*reshape( cb-mm, len, 1 );
        wx             =   repmat( vec, 1, cls_num );
        dis            =   sum( (wx - codewords).^2 );        
        [md, idx]      =   min(dis);
        PCA_idx(cnt)   =   idx;
        cnt            =   cnt + 1;
    end
    
end
        