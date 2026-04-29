function [HurstXY,debugob] = bbctv(Im,          ...
    				   firs,        ...
    		                   rr,          ...
    				   ds,          ...
   			           weight_flag, ...
                                   con_flag,    ...
                                   lam,         ...
                                   D,           ...
    				   max_iter)

% ----------------------------------------------------------------------
% function [HurstXY,debugob] = bbctv(Im,          ...
%    				     firs,        ...
%    		                     rr,          ...
%    				     ds,          ...
%   			             weight_flag, ...
%                                    con_flag,    ...
%                                    lam,         ...
%                                    D,           ...
%    				     max_iter)
% ----------------------------------------------------------------------
% Use tv regularistaion to estimate Hurst exponent via LS
% ----------------------------------------------------------------------
% Note downsampling (resize energy at all levels to 
% [M/ds, N/ds]) rather than [M, N]
% ----------------------------------------------------------------------
% . Im = 2-d array/image
% . firs.biort, firs.qshift = filter options for dtcwt function
% . rr(1) = coarsest scale cut-off (typically = 2)
% . rr(2) = finest scale cut-off (typically = N_levels-1)
% . ds = downsampling factor
% . weight_flag = perform weighted least squares (nb this is different to 
% . ols only when ds not = 1)
% . con_flag = perform box constrained version estimation
% . lam = tv regularisor parameter
% . D = penalty matrix (use D = penalty2(M, N, ds) for 2d differencer); 
% . maxiter = maximum iterations used  
% ----------------------------------------------------------------------
% . HurstXY = output local Hurst estimate
% . debugob = output debug object (not used at the moment)


%% Preliminaries -------------------------------------------------------
biort = firs.biort;
qshift = firs.qshift;

meanIm   = mean(Im(:));
Im       = Im-meanIm;
[M, N]   = size(Im);
N_levels = log2(min(M,N));
p        = 1; 

%% DTCWT ---------------------------------------------------------------
[Yl, Yh, Yscale] = dtwavexfm2(Im, N_levels, biort, qshift);

%% Compute energies ----------------------------------------------------
for d = 1:6
    for r=1:N_levels
% --- downsampling 
        ye{d}{r} = imresize(abs(Yh{r}(:,:,d)).^p, ...
	                    [M/ds,N/ds],          ...
			    'bilinear');
    end
end

Exy_t = zeros(M*N/(ds*ds), 6, N_levels);
C = M*N/(ds*ds);
for r = 1:N_levels
    c = 0;
    for th=1:6
        sY = ye{th}{r};
        tmp = log2(sY.^(3-p)+1e-12);
        Exy_t(c+1:c+C, th, N_levels-r+1) = tmp(:);
    end
end

%% Log least squares fit -----------------------------------------------
rr = rr(1):rr(2);
R(:,1)  = ones(length(rr), 1);
R(:,2)  = rr';

%% set up constraints (for beta) ---------------------------------------
% === x_lo and x_hi are the constraints
% --- first, set up generous constraints; equiv to +/- inf 
x_lo = -1000000*ones(2*M*N/(ds*ds), 1);
x_hi =  1000000*ones(2*M*N/(ds*ds), 1);
% --- now only constrain the slope to be in [-4,-2]; 
% (this is equiv to constraining H in [0,1]):
if con_flag
  x_lo(2:2:end) = -4;
  x_hi(2:2:end) = -2;
end

%% estimation ----------------------------------------------------------
% --- nb cols of A (assuming R always has 2 columns i.e. that there are two
% parameters per pixel, namely an intercept and a slope) ---
nbcols_A = 2*M*N/(ds*ds);  
dx = zeros(max_iter,th);

% --- block matrix version of R ----------------------------------------
A = kron(speye(M*N/(ds*ds)), R); 

xbar = 1/2 * (x_lo + x_hi);
[m,n] = size(A);
for th = 1:6
% --- initialise constraints parameters --------------------------------
  stepeps = ones(n,1);
  count   = 0;
  invC    = diag(4./(x_hi - x_lo).^2);
  mask    = ones(n,1);

% --- stack energy into vector -----------------------------------------
  y = squeeze(Exy_t(:,th,rr))';
  y = y(:);
    
% --- initial estimate -------------------------------------------------
  x = inv(A'*A + lam*D'*D + invC) * (A'*y + invC*xbar);
% --- start iterations -------------------------------------------------
  while sum(mask) && (count < max_iter)
    count = count+1;
% --- tv reg update ----------------------------------------------------
    Om = diag(1./(sqrt(abs(D*x))+1e-5));
% --- solution update --------------------------------------------------
    x = inv(A'*A + lam*D'*Om'*Om*D + invC) * (A'*y + invC*xbar);
% --- constraint updates -----------------------------------------------
    mask    = (x >= x_hi) | (x <= x_lo);
    stepeps = 1/(1+count/10) * stepeps;
    invC    = diag(~mask + mask./stepeps) .* invC;
% --- diagnostics ------------------------------------------------------
    x_tmp = x;
    dx(count,th) = norm(x-x_tmp);
  end
  Hurstxy_t(th,:) = -1 - x(2:2:end)'/2;
  HurstXY_t(:,:,th) = reshape(Hurstxy_t(th,:),M/ds,N/ds);
end

Hurstxy = mean(Hurstxy_t,1);
HurstXY = reshape(Hurstxy,M/ds,N/ds);

% --- "upsample" back to original size ---
HurstXY = imresize(HurstXY, [M,N], 'bilinear');

debugob.dx = dx;
return
%% ---------------------------------------------------------------------



