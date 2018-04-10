function[x,w] = hermite_gauss_quadrature(N)
% Returns the N-point hermite-gauss quadrature rule over the real line, integrating with respect to weight function
%
%  w(x) = exp(-x^2)

a = zeros([N 1]);
b = ones([N 1]);
b(1) = gamma(1/2);
b(2:end) = (1:(N-1))./2;

J = spdiags([sqrt([b(2:N);0]) a(1:N) sqrt(b(1:N))], -1:1, N, N);
[V, x] = eig(full(J), 'vector');

w = (V(1,:).^2).';
w = w*sqrt(pi)/sum(w);
