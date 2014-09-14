function d_im = IDR_Shrinkage( d_im, hist_refh, hist_refv, par)

[m, n] = size(d_im);
trans = @(X) 1/sqrt(m*n)*fft2(X);
itrans = @(X) sqrt(m*n)*ifft2(X);
eigDtD = abs(psf2otf([1,-1],[m,n])).^2 + abs(psf2otf([1;-1],[m,n])).^2;

miu = par.miu; %5e-2;
d_im0 = d_im;
for iter = 1:10
    % Update dh and dv
    [dh, dv] = Ltrans2(d_im);
    sign_dh = sign(dh);
    sign_dv = sign(dv);
    dh0 = 255*histeq(abs(dh)/255,hist_refh);
%     dh0 = min(abs(dh), dh0);
    dh  = sign_dh .* dh0;
    dv0 = 255*histeq(abs(dv)/255,hist_refv);
%     dv0 = min(abs(dv), dv0);
    dv  = sign_dv .* dv0;
    
    % Update d_im
    x_d = miu * d_im0 + Lforward2(dh, dv);
    d_im = real( itrans( trans(x_d) ./ (miu + eigDtD) ) );
    d_im = max(d_im, 0);
    d_im = min(d_im, 255);
end

