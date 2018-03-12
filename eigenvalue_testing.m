clear
close all

a = 1; b = 0; c = 1; d = 0;

A = [a b; c d];

%A = randn(2);

N = 100;

[lambda, k] = positive_eigenvalues(A, N);

phi =  @(x,n) k(n,1)*cos(x*sqrt(lambda(n))) + k(n,2)*sin(x*sqrt(lambda(n)));
phid = @(x,n) -k(n,1)*sqrt(lambda(n))*sin(x*sqrt(lambda(n))) + k(n,2)*sqrt(lambda(n))*cos(x*sqrt(lambda(n)));

% Test boundary conditions
phi0 = zeros([N 1]);
phi1 = zeros([N 1]);
for n = 1:N
  phi0(n) = A(1,1) * phi(0, n) + A(1,2) * phid(0,n);
  phi1(n) = A(2,1) * phi(1, n) + A(2,2) * phid(1,n);
end

% Test orthogonality
[x,w] = gauss_quadrature(2*N, 0, 1);

V = zeros([length(x) N]);
for n = 1:N
  V(:,n) = sqrt(w).*phi(x,n);
end
G = V'*V;

% M = 6;
% set(plot(x, diag(1./sqrt(w))*V(:,1:M)), 'linewidth', 4);
