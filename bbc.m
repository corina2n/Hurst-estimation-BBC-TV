function [HurstXY, dob] = bbc(Im,          ...
                              firs,        ...
    		              rr,          ...
    			      ds,          ...
   			      weight_flag, ...
                              con_flag)

% ----------------------------------------------------------------------
% function [HurstXY, debugob] = dtcwtConHurst(Im,          ...
%                                             firs,        ...
%     					      rr,          ...
%    					      ds,          ...
%   					      weight_flag, ...
%                                             con_flag)
% ----------------------------------------------------------------------
% Use dtcwt to estimate Holder exponent via OLS with option (via
% con_flag) to do constrained ols and option (via weight_flag) to do
% weighted least squares to take into account different numbers of
% samples at each scale level (when ds not= 1)
% ----------------------------------------------------------------------
% . Im = 2-d array/image
% . firs.biort, firs.qshift = filter options for dtcwt function
% . rr(1) = coarsest scale cut-off (typically = 2)
% . rr(2) = finest scale cut-off (typically = N_levels-1)
% . ds = downsampling factor
% . weight_flag = perform weighted least squares (nb this is different to 
% . ols only when ds not = 1)
% . con_flag = perform box constrained version estimation
% ----------------------------------------------------------------------
% . HurstXY = output local Hurst estimate
% . debugob = output debug object (not used at the moment)
% ----------------------------------------------------------------------
% . Note the downsampling (we resize energy at all levels to [M/ds,
% N/ds] rather than [M,N] 
% ----------------------------------------------------------------------
% . Note the two equivalent ways of setting up the OLS solution.  The
% 1st is the usual way; the 2nd is by introducing a huge matrix X s.t. y
% = Ax where   y is the energy of the wavelet coeffs and where the
% vector x  contains the estimated slope of the wavelet spectra (in its
% even indexed elements). (Its odd indexed elements are the
% intercepts.) 
%
% *** Caveat: when ds>1 this will use bilinear interp to rescale
% estimated Hurst back at the end (depending on whether you comment out
% the line following:
% --- "upsample" back to original size ---

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
    for r=1:log2(ds)
       ss = max(1, ds*2.^(N_levels-r)./N);
       ye{d}{r} = blkproc(abs(Yh{r}(:,:,d)).^p, [ss, ss], 'mean2');
    end

    for r=log2(ds)+1:N_levels
% --- downsampling (if ds not = 1)
        ye{d}{r} = imresize(abs(Yh{r}(:,:,d)).^p, ...
	                    [M/ds, N/ds],         ...
			    'bilinear');
    end
end

Exy_t = zeros(M*N/(ds*ds), 6, N_levels);
C = M*N/(ds*ds);
for r = 1:N_levels
    c = 0;
    for th=1:6
        sY = ye{th}{r};
        %tmp = log2(sY.^2+1e-12);
        %tmp = log2(sY+1e-12);
        tmp = log2(sY.^(3-p)+1e-12);
        Exy_t(c+1:c+C, th, N_levels-r+1) = tmp(:);
    end
end

%% Log least squares fit -----------------------------------------------
rr = rr(1):rr(2);
R(:,1)  = ones(length(rr), 1);
R(:,2)  = rr';

if weight_flag==1
  w = 1./(2.^(-rr)); 
  W0 = diag(w);
  %w = max(1, ds*2.^(0:log2(N)-1)./N);
  %w = w(rr);
  %w(end) = w(end-1);
  %W0 = diag(w);
% --- experimental -----------------------------------------------------
  elseif weight_flag==2
    w = max(1, ds*2.^(0:log2(N)-1)./N);
    w = w(rr);
    W0 = diag(w);
  else
% --- unweighted case
    w = ones(length(rr),1);
    W0 = eye(length(rr));
end
W = kron(speye(M*N/(ds*ds)), W0); 

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
Hurstxyt = zeros(6, M*N/(ds*ds));
maxiter = 100;
verb = 0;
for th = 1:6
%  y0 = squeeze(Exy_t(:, th, rr))';
  y = squeeze(Exy_t(:, th, rr))';
% --- stack into vector 
%  y = y0(:); 
  y = y(:); 
% --- block matrix version of R
  A = kron(speye(M*N/(ds*ds)), R); 
% --- if constriant flag is on, use constrained ols; otherwise use ols:
  if con_flag
    x = constr_ols(A, y, x_lo, x_hi, W, maxiter, verb); 
  else
%    x2 = A \ y; 
% --- weighted least squares solver
    x = lscov(A, y, diag(W));
%    x3 = inv(R'*R)*R'*y0;
  end
% --- rearrange beta estimate vector into Hurst estimate vector
  Hurstxy_t(th,:) = -1 - x(2:2:end)'/2;
end

% --- take mean over orientations
Hurstxy = mean(Hurstxy_t, 1);
% --- reshape Hurst estimates into image 
HurstXY = reshape(Hurstxy, M/ds, N/ds);

% --- "upsample" back to original size ---
HurstXY = imresize(HurstXY, [M,N], 'bilinear');
% --- debugging object:
dob = [];
return

%dob.W = W;
%dob.A = A;
%dob.R = R;
%dob.y = y;
%dob.y0 = y0;
%dob.E = Exy_t;
%dob.ye = ye;
%dob.rr = rr;
%dob.x2 = x2;
%dob.x3 = x3;
%dob.x = x;

%% end of file ---------------------------------------------------------


