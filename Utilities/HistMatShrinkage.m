function [dh, dv] = HistMatShrinkage( d_im, hist_refh, hist_refv)

[dh, dv] = Ltrans2(d_im);
sign_dh = sign(dh);
sign_dv = sign(dv);
dh0 = 255*histeq(abs(dh)/255,hist_refh);
%     dh0 = min(abs(dh), dh0);
dh  = sign_dh .* dh0;
dv0 = 255*histeq(abs(dv)/255,hist_refv);
%     dv0 = min(abs(dv), dv0);
dv  = sign_dv .* dv0;
