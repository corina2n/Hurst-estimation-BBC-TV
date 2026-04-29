fprintf('If you use this code please cite:\n ');
fprintf('J. D. B. Nelson; C. Nafornita; A. Isar,\n');
fprintf('\"Semi-local scaling exponent estimation with box-penalty constraints and total-variation regularisation,\" \n');
fprintf('in IEEE Transactions on Image Processing, doi: 10.1109/TIP.2016.2551365,\n');
fprintf('URL: http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=7448417 \n');


%% === comments --------------------------------------------------------
% [1] generate fractal Brownian surface with piecewise constant varying
%Hurst parameter 
%
% [2] estimate Hurst using DTCWT method with configurable combination
% of: energy downsampling; semi-local weights; bound constraints on
% beta; total variation smoothing

%% === preliminaries ---------------------------------------------------
close all
clear all

%% === make sure to tell Matlab where the dtcwt toolbox is -------------
% !!! CAVEAT: this assumes presence of INRIA's fraclab and Nick
% Kingsbury's DTCWT toolbox 
addpath('dtcwt', genpath('fraclab'))

%% === size parameters -------------------------------------------------
N = 256;             
% --- height -----------------------------------------------------------
M = N;              
N_levels = log2(min(M,N));

%% === DTCWT parameters ------------------------------------------------
% --- level 1 FIR ------------------------------------------------------
firs.biort  = 'near_sym_b_bp'; 
% --- Q-shift FIR ------------------------------------------------------
firs.qshift = 'qshift_b_bp';   
% --- coarsest level ---------------------------------------------------
j0 = 2;
% --- finest level -----------------------------------------------------
j1 = N_levels-2;
J = [j0, j1];

%% === experimental parameters -----------------------------------------
% --- true Hurst -------------------------------------------------------
H = 0.3;

%% === methodological parameters ---------------------------------------
% --- wavelet energy downsampling --------------------------------------
ds = 16;
% --- use weights ------------------------------------------------------
wflag = 1;
% --- use constraints --------------------------------------------------
cflag = 1;

%% === simulate fBf ----------------------------------------------------
%y = fastfBm2D([M N], H);
fprintf('\n------------------------------------\n')
fprintf('synthesising fBf... ')
[y,Hxy] = synthmbf([M, N]);
fprintf('done!\n')
% --- ols --------------------------------------------------------------

fprintf('estimating Hurst with ols... ')
[Hols, dob] = bbc(y, firs, J, ds, 0, 0);
fprintf('done!\n')

% --- BBC-WLS ----------------------------------------------------------
fprintf('estimating Hurst with bbc... ')
[Hbbcwls, dob] = bbc(y, firs, J, ds, wflag, cflag);
fprintf('done!\n')
	
% --- total variation --------------------------------------------------
D = penalty2(M, N, ds); 
maxiter = 15;
lam = 0.5;

fprintf('estimating Hurst with bbctv... ')
[Hbbctv, dob] = bbctv(y, firs, J, ds, wflag, cflag, lam, D, maxiter);
fprintf('done!\n')

subplot(221)
imagesc(y)
%title(fprintf('fBf, Hurst=%.2f', H))
title('fBf')

subplot(222)
imagesc(Hols, [0 1])
title('OLS')

subplot(223)
imagesc(Hbbcwls, [0 1])
title('BBC-WLS')

subplot(224)
imagesc(Hbbctv, [0 1])
title('BBC-TV')

%mean2(Hols)
%mean2(Hbbcwls)
%mean2(Hbbctv)

fprintf('ols error:    %.4f\n', mean2(abs(Hxy-Hols)))
fprintf('bbcwls error: %.4f\n', mean2(abs(Hxy-Hbbcwls)))
fprintf('bbctv error:  %.4f\n', mean2(abs(Hxy-Hbbctv)))




