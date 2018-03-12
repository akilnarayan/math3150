function[lambda, k] = zero_eigenvalues(A)
% zero_eigenvalues -- Determines if 0 is an eigenvalue for an ODE
%
% [lambda, k] = zero_eigenvalues(A)
%
%   Determines if the value lambda = 0 is an eigenvalue for the boundary value
%   problem
%
%          y''(x) + lambda*y(x) = 0,   0 < x < 1
%
%           A(1,1)*y(0) + A(1,2)*y'(0) = 0,
%           A(2,1)*y(1) + A(2,2)*y'(1) = 0.
%
%   If so, then the output lambda is set to 0, and k is a 1 x 2 vector so that 
%
%     y(x) = k(1) + k(2) * x
%
%   is the corresponding orthonormal eigenfunction.

lambda = [];
k = zeros([0 2]);

[a,b,c,d] = deal( A(1,1), A(1,2), A(2,1), A(2,2) );

if ( a == 0 ) && ( c == 0 )
  lambda = 0;
  k = [1 0];
else %if ( a ~=0 ) && ( c ~= 0) 

  if abs( b/a - d/c - 1 ) < 1e-13

    k = [A(1,2)  -A(1,1)];

    knorm = (k(2) - k(1)/2)^2 + k(1)^2 / 12;

    k = k/sqrt(knorm);

  end

end
