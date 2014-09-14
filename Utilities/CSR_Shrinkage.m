function  [im_out]  =  CSR_Shrinkage( im, par, alpha, beta, Tau1, Dict, flag )
[h w ch]   =   size(im);
A          =   Dict.D0;
PCA_idx    =   Dict.cls_idx;
s_idx      =   Dict.s_idx;
seg        =   Dict.seg;
b          =   par.win;
b2         =   b*b*ch;
Y          =   zeros(b2, size(alpha,2), 'single' );

idx          =   s_idx(seg(1)+1:seg(2));
tau1         =   par.tau1;
if  flag==1
    tau1    =   Tau1(:, idx);
end
Y(:, idx)    =   A'*(soft( alpha(:,idx)-beta(:,idx), tau1 ) + beta(:,idx));

for   i  = 2:length(seg)-1   
    idx    =   s_idx(seg(i)+1:seg(i+1));    
    cls    =   PCA_idx(idx(1));
    P      =   reshape(Dict.PCA_D(:, cls), b2, b2);
    tau1         =   par.tau1;
    if  flag==1 
        tau1    =   Tau1(:, idx);
    end
    Y(:, idx)    =   P'*(soft( alpha(:,idx)-beta(:,idx), tau1 ) + beta(:,idx));
end

s     =  par.step;
N     =  h-b+1;
M     =  w-b+1;
r     =  [1:s:N];
r     =  [r r(end)+1:N];
c     =  [1:s:M];
c     =  [c c(end)+1:M];
N     =   length(r);
M     =   length(c);
im_out   =  zeros(h,w);
im_wei   =  zeros(h,w);
k        =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        im_out(r-1+i,c-1+j)  =  im_out(r-1+i,c-1+j) + reshape( Y(k,:)', [N M]);
        im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + 1;
    end
end

im_out  =  im_out./(im_wei+eps);
return;



function Y   =  CSR_threshold( x, mu, tau1, tau2 )
m1   =  (mu>=0).*mu;
m2   =  (mu<0).*mu;

y1   =  ( x < (-tau1-tau2) ).*( x+tau1+tau2 );
y3   =  ( x>(tau1-tau2) & x<(tau1-tau2+m1) ).*( x-tau1+tau2 );
y4   =  ( x>=(tau1-tau2+m1) & x<=(tau1+tau2+m1) ).*m1;
y5   =  ( x>(tau1+tau2+m1) ).*( x-tau1-tau2 );
Y1   =  (y1 + y3 + y4 + y5).*(mu>=0);

y1   =  ( x<(m2-tau1-tau2) ).*( x+tau1+tau2 );
y2   =  ( x>=(m2-tau1-tau2) & x<=(tau2+m2-tau1) ).*m2;
y3   =  ( x>(tau2+m2-tau1) & x<(tau2-tau1) ).*( x+tau1-tau2 );
y5   =  ( x>(tau1+tau2) ).*( x-tau1-tau2 );
Y2   =  (y1 + y2 + y3 + y5).*(mu<0);

Y    =  Y1 + Y2;
return;


function y   =  CSR_soft(x, mu, tau1, tau2)
tau1    =   2*tau1;
x       =   x./(1+tau1) - mu;
tau     =   tau2./(1+tau1);
y       =   sign( x ).*max(abs(x)-tau,0) + mu;
return;
