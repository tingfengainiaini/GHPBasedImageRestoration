function  Y  =  FWD_CPCA( im, CPCA_D, mX, par )

b    =  par.win;
st   =  par.step;

[h w ch] =  size(im);

num     =  ceil(h/st)*ceil(w/st);
Y       =  zeros(b*b*ch,num);
cnt     =  1;
for  row  =  1:st:h

    if  (row+b-1) > h
        row   =  h-b+1;
    end       
    
    for  col  =  1:st:w
        
        if  (col+b-1) > w
            col   =  w-b+1;
        end                             
        
        cb      =  im(row:row+b-1, col:col+b-1, :);
        cb(:,:,1)   =  cb(:,:,1)';
        cb(:,:,2)   =  cb(:,:,2)';
        cb(:,:,3)   =  cb(:,:,3)';
        
        P           =  reshape( double(CPCA_D(:, cnt)) , (b^2)*ch, (b^2)*ch );
        
        Y(:,cnt)    =  P*(cb(:) - double(mX(:, cnt)));
        
        cnt   =  cnt + 1;
        
    end
end
        
