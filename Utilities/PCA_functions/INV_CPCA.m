function  im  =  INV_CPCA( Y, CPCA_D, mX, par )

b    =  par.win;
st   =  par.step;
h    =  par.h;
w    =  par.w;
ch   =  3;
im   =  zeros(h, w, ch);
wei  =  zeros(h, w);

cnt     =  1;
for  row  =  1:st:h

    if  (row+b-1) > h
        row   =  h-b+1;
    end       
    
    for  col  =  1:st:w
        
        if  (col+b-1) > w
            col   =  w-b+1;
        end                             
       
        
        P     =  reshape( double(CPCA_D(:, cnt)), (b^2)*ch, (b^2)*ch);
        cb    =  P'*Y(:, cnt) + double(mX(:, cnt));
       
        cb    =  reshape(cb, b, b, ch);
        cb(:,:,1)  =  cb(:,:,1)';
        cb(:,:,2)  =  cb(:,:,2)';
        cb(:,:,3)  =  cb(:,:,3)';
        
        im(row:row+b-1, col:col+b-1, :)    =  im(row:row+b-1, col:col+b-1, :) + cb;
        wei(row:row+b-1, col:col+b-1)      =  wei(row:row+b-1, col:col+b-1) + 1;
        
        
        cnt   =  cnt + 1;
        
    end
end

im(:,:,1)    =  im(:,:,1)./wei;
im(:,:,2)    =  im(:,:,2)./wei;
im(:,:,3)    =  im(:,:,3)./wei;