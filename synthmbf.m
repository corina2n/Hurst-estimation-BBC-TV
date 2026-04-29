function [B,Hxy] = synthmbf(sz)
% Synthesises very simple Multifractional Brownian Surface

%% === quantise Hxy into 4-many levels ---------------------------------
M = sz(1);
N = sz(2); 
B     = zeros(M,N);
qbin  = false(M,N,4);
H = [0.2, 0.4, 0.6, 0.8];

% --- partition support ------------------------------------------------
Hxy = zeros(M,N);
Hxy(1:M/2, 1:N/2)     = H(1);
Hxy(M/2+1:M, 1:N/2)   = H(2);
Hxy(1:M/2, N/2+1:N)   = H(3);
Hxy(M/2+1:M, N/2+1:N) = H(4);

% --- variance of Wx, Wy equal to 1 ------------------------------------
W   = fft2(randn(2*M,2*N));
Wx  = fft(randn(2*M,1))/sqrt(2*M); 
Wy  = fft(randn(1,2*N)).'/sqrt(2*N);

%% --- synthesise and stitch together ----------------------------------
for q = 1:4
% --- simulate fBs with H = q/Q ----------------------------------------
    Imq = fastfBm2D([M N], H(q), W, Wx, Wy);
% --- find pixels of Hxy where Holder = q/Q... -------------------------
    qbin(:,:,q) = Hxy==H(q);
% --- set those pixels of B to Imq -------------------------------------
    B(qbin(:,:,q)) = Imq(qbin(:,:,q));
end











