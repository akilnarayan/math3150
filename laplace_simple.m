% Computes solution to 1D Laplace's equation with no forcing...i.e., super easy.

clear
close all

N = 128; % Number of eigenvalue/eigenfunction pairs to compute
L = 1; % Length of the domain
k = 0.01; % u_t = k u_xx

% Initial data -- should be a function defined in a vectorized way
f = @(x) x.*(1-x).*exp(cos(4*abs(x-0.5)));

% Boundary conditions
%  a * u(0) + b * u'(0) = v(1)
%  c * u(L) + d * u'(L) = v(2)
a = 1; b = 0;
c = 1; d = 0;
v = [0; 3];

tol = 1e-12;

% Preprocessing -- normalize boundary data
temp = norm([a b]); assert (temp > tol, 'Must specify nontrivial boundary conditions');
a = a/temp; b = b/temp; v(1) = v(1)/temp;
temp = norm([c d]); assert (temp > tol, 'Must specify nontrivial boundary conditions');
c = c/temp; d = d/temp; v(2) = v(2)/temp;

Abc = [b a; c*L+d c];

% Compute equilibrium solution to handle non-homogenous boundary conditions
c = pinv(Abc)*v;
if norm(Abc*c - v) > 1e-12
  error('Boundary conditions do not satisfy solvability condition');
else
  ff = @(x) f(x) - c(1)*x - c(2);
end

M = 2*N; % Number of quadrature points for approximation of integrals
[x,w] = gauss_quadrature(M, 0, L);

%%%%%% Visualization options
Nviz = 8;
lineprops = {'linewidth', 3};
labelprops = {'fontsize', 16, 'fontweight', 'b', 'interpreter', 'latex'};
axesprops = {'fontsize', 16, 'fontweight', 'b'};
%%%%%%

figure;
set(plot(x, c(1)*x + c(2), 'r'), lineprops{:});
set(xlabel('$\mathbf{x}$'), labelprops{:});
set(ylabel('$\mathbf{u(x)}$'), labelprops{:});
set(gca, axesprops{:});
