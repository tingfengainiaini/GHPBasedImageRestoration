%----------------------------------------------------
% Blur: 9x9 Uniform kernel; AGWN std Varance = 2^0.5
% Data: May 20th, 2010
% Author: Weisheng Dong, wsdong@mail.xidian.edu.cn
%----------------------------------------------------
function  Image_Denoising(method, nSig, Out_dir, In_dir )
par.method    =   method;
pre           =   'NCSR_';
par.nSig      =   nSig;
par.iters     =   1;
par.eps       =   1e-8;
par.method    =   method;
par.cls_num   =   70;   % 70

if nSig<=30
    par.c1        =   0.64;         % 0.57
else    
    par.c1        =   0.64;     
end
par.lamada    =   0.23;  % 0.36
par.hh        =   12;

if nSig <= 10
    par.c1        =   0.56;   
    par.lamada    =   0.22;
    
    par.sigma     =   1.7;  
    par.win       =   6;
    par.nblk      =   13;           
    par.hp        =   75;
    par.K         =   3;
    par.m_num     =   240;
elseif nSig<=15
    par.c1        =   0.59;   
    par.lamada    =   0.22;
    
    par.sigma     =   1.8;
    par.win       =   6;
    par.nblk      =   13;           
    par.hp        =   75;
    par.K         =   3;
    par.m_num     =   240;    
elseif nSig <=30
    par.sigma     =   2.0;
    par.win       =   7;
    par.nblk      =   16;
    par.hp        =   80;       %  80
    par.K         =   3;
    par.m_num     =   250;
elseif nSig<=50
    par.c1        =   0.64;    
    
    par.sigma     =   2.4;
    par.win       =   9;
    par.nblk      =   18;
    par.hp        =   90;
    par.K         =   4;
    par.m_num     =   300;
    par.lamada    =   0.26;    
else
    par.sigma     =   2.4;
    par.win       =   8;
    par.nblk      =   20;
    par.hp        =   95;
    par.K         =   4;
    par.m_num     =   300;
    par.lamada    =   0.26;
end

fpath         =   fullfile(In_dir, '*.png');
im_dir        =   dir(fpath);
im_num        =   length(im_dir);
cnt           =   0;
sum_psnr      =   0;
sum_ssim      =   0;

time0         =   clock;
fn_txt        =   strcat( pre, 'PSNR_SSIM.txt' ); 
fd_txt        =   fopen( fullfile(Out_dir, fn_txt), 'wt');

% par.maskstr     =  'nature_segm.png'; % 'orig_10_6_10gseg_segm.png'; % 'Monarch_segm.png'; % 'house_segm.png';
% par.maskstr     =  'orig_10_6_10gseg_segm.png'; % ''; % 'Monarch_segm.png'; % 'house_segm.png';
par.maskstr         =   'tmp_segm.png';
global      SHOW

SHOW        =   1;

for i = 2:2%im_num 
    
%     par.I        =   double( imread(fullfile(In_dir, im_dir(i).name)) );
%     par.I        =   double( rgb2gray( imread(fullfile(In_dir, im_dir(i).name)) ) );
%     par.I        =   par.I; %(1:128,128+1:128+128);%(128+1:256, 64+1:256);
%     Noise        =   nSig*Gen_noise(In_dir, im_dir, i);
%     par.nim      =   par.I + Noise;   % randn(size(par.I));
    load('house.png.mat');
    par.I   =   I;
    par.nim     =   nim;
    [par.hist_dh, par.hist_dv]   =   HistEst(par.nim, nSig);
%     [par.hist_dh, par.hist_dv, par.par_h, par.par_v, par.mask] = HMSeg(par.nim, par.maskstr, nSig);
    
    [im PSNR SSIM]   =   Centralized_SR_Denoising_IDR2( par );      %_IDR2 
    sum_psnr    =  sum_psnr + PSNR;
    sum_ssim    =  sum_ssim + SSIM;    
    
    fname            =   strcat(pre, im_dir(i).name);    
    imwrite(im./255, fullfile(Out_dir, fname));
    disp( sprintf('%s: PSNR = %3.2f  SSIM = %f\n', im_dir(i).name, PSNR, SSIM) );
    fprintf(fd_txt, '%s :  PSNR = %2.2f  SSIM = %2.4f\n', im_dir(i).name, PSNR, SSIM);
    cnt   =  cnt + 1;
end
fprintf(fd_txt, '\n\nAverage :  PSNR = %2.2f  SSIM = %2.4f\n', sum_psnr/cnt, sum_ssim/cnt);
fclose(fd_txt);
disp(sprintf('Total elapsed time = %f min\n', (etime(clock,time0)/60) ));
return;


function nim  =  Gen_noise( In_dir, im_dir, i )
randn('seed',0);
for ii=1:i
    im        =   rgb2gray( imread(fullfile(In_dir, im_dir(ii).name)) );
    nim       =   randn(size(im));
end
return;

