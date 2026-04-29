function x = constr_ols(A, b, x_lo, x_hi, W, maxiter, verb)
% 
% solve ||Ax - b||_2 s.t. x_lo <= x <= x_hi
%
% using method on p4 of Mead and Renaut
% Algorithm 1 Solve least squares problem with box constraints as quadratic constraints
%
% J. D. B. Nelson (UCL) and C. Nafornita (Politehnica University
% Timisoara); June 2015

%% === preliminaries ---------------------------------------------------
[M,N] = size(A);
% --- centre of feasible region
xbar = 1/2 * (x_lo + x_hi);
if verb
  fprintf('\n ----------------------------------------------------\n') 
end

%% === initialization --------------------------------------------------
stepeps = ones(N,1);
count   = 0;
invC    = diag(4./(x_hi - x_lo).^2);
mask    = ones(N,1);

%% === iterations ------------------------------------------------------
while sum(mask) && (count < maxiter)
  count   = count+1;
  x       = inv(A'*W*A + invC)*(A'*W*b + invC*xbar);
  mask    = (x >= x_hi) | (x <= x_lo);
  stepeps = 1/(1+count/10) * stepeps;
  invC    = diag(~mask + mask./stepeps) .* invC;

  if verb 
    disp(x'), disp(mask'), disp(diag(invC)'), fprintf('\n')
  end
end


 
