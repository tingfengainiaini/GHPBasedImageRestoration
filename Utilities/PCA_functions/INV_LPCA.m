function  im  =  INV_LPCA( Y, PCA_idx, PCA_D, par )

b    =  par.win;
st   =  par.step;
h    =  par.h;
w    =  par.w;

im   =  zeros(h, w);
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
               
        
        P       =  reshape(PCA_D(:, PCA_idx(cnt)), b^2, b^2);
        cb      =  P'*Y(:, cnt);
        
        
        im(row:row+b-1, col:col+b-1)    =  im(row:row+b-1, col:col+b-1) + reshape(cb, b, b)';
        wei(row:row+b-1, col:col+b-1)   =  wei(row:row+b-1, col:col+b-1) + 1;
        
        
        cnt   =  cnt + 1;
        
    end
end

im    =  im./wei;