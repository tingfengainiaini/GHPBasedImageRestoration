function  Y  =  FWD_PCA( im, PCA_D, mX, par )

b    =  par.win;
st   =  par.step;

[h  w]  =  size(im);

num     =  floor(h*w/st);
Y       =  zeros(b*b,num);
cnt     =  1;
for  row  =  1:st:h

    if  (row+b-1) > h
        row   =  h-b+1;
    end       
    
    for  col  =  1:st:w
        
        if  (col+b-1) > w
            col   =  w-b+1;
        end                             
        
        cb      =  im(row:row+b-1, col:col+b-1)';
        
        P       =  reshape( double(PCA_D(:, cnt)) , b^2, b^2);
        
        Y(:,cnt)   =  P*(cb(:) - double(mX(:, cnt)));
        
        cnt   =  cnt + 1;
        
    end
end
        
