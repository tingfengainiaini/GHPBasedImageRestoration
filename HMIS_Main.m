%----------------------------------------------------
% Texture Enhanced Image Denoising via Gradient Histogram Preservation
% Author: Wangmeng Zuo, Lei Zhang, Chunwei Song, David Zhang
%----------------------------------------------------

% Please inpute addpath( genpath('Utilities') ) in command window first.
% I is the noiseless image, which is 256-level gray image
% nSig is the noise standard deviation.
% this function here is just a interface for testing.
function    dim  = HMIS_Main( I, nSig)

nim    =   Add_noise( I, nSig );

dim    =    Image_Denoising_HMIS( nSig, nim, I );  

imwrite(dim/255, 'Results\denoise.png');
imwrite(I, 'Results\noiseless.png');
imwrite(nim/255, 'Results\noise.png');

 