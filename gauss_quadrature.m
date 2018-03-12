function[x,w] = gauss_quadrature(N, left, right)
% (Legendre) Gaussian quadrature
%
% [x,w] = gauss_quadrature(N)
%
%   Computes the N-point Gaussian quadrature rule associated to the integral
%
%    \int_{-1}^1 f(x) dx = \sum_{n=1}^N w(n) f(x(n))
%
% [x,w] = guass_quadrature(N, left, right)
%
%   Computes the Gauss quadrature rule for integration over [left, right].

assert(N >= 1);

% Recurrence coefficients: Special case of Jacobi polynomials with alpha = beta = 0
a = zeros([N 1]);
b = ones([N 1]);

b(1) = 1;
if N > 1
  ns = 1:(N-1);
  b(2:end) = ns.^2./(4*ns.^2 - 1);
end

J = spdiags([sqrt([b(2:N);0]) a(1:N) sqrt(b(1:N))], -1:1, N, N);
[V, x] = eig(full(J), 'vector');

w = (V(1,:).^2).';
w = w*2/sum(w);

if nargin == 3
  assert(left < right);
  x = (x + 1)/2 * (right - left) + left;
  w = w/2*(right-left);
end
