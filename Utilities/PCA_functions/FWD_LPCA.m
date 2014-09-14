function  Y  =  FWD_LPCA( im, PCA_idx, PCA_D, par )

b    =  par.win;
st   =  par.step;

[h  w]  =  size(im);

Y       =  zeros(b*b, floor(h*w/st));
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
        
        P       =  reshape(PCA_D(:, PCA_idx(cnt)), b^2, b^2);
        
        Y(:,cnt)   =  P*cb(:);
%         Y          =  [Y  P*cb(:)];
        
        cnt   =  cnt + 1;
        
    end
end
        
