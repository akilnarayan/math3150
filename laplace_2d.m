% Computes an approximate solution to the one-dimensional Laplace equation
% (ODE) on [-1,1]^2 with periodic boundary conditions:
%
%    u_xx + u_yy = f

clear
close all

N = 1000;                % spatial points in 1D

% Right-hand side
f = @(x,y) exp(cos(35*x)).*sin(sin(4*y)) + exp(cos(5*pi*y)).*exp(sin(pi*x)) - cosh(2) - cosh(1);

x = linspace(-1, 1, N+1).';
x(end) = [];
h = 2/N;

[xx,yy] = ndgrid(x, x);

frhs = f(xx,yy);

% Form 1D differentiation operator
A = zeros(N);
A = A - 2/h^2 * diag(ones([N 1])) + 1/h^2 * diag(ones([N-1 1]), 1) + 1/h^2 * diag(ones([N-1 1]), -1);
A(1,N) = 1/h^2;
A(N,1) = 1/h^2;

% Resolve kernel
A = A + 1/N;

u = sylvester(A, A, frhs);

labelprops = {'fontsize', 16, 'fontweight', 'b', 'interpreter', 'latex'};

figure;
subplot(1,2,1);
surf(xx, yy, frhs);
shading interp
set(title('$\mathbf{f(x)}$'), labelprops{:});
set(gca, 'fontsize', 20, 'fontweight', 'b');

subplot(1,2,2);
surf(xx, yy, u);
shading interp
set(title('$\mathbf{u(x,t)}$'), labelprops{:});
set(gca, 'fontsize', 20, 'fontweight', 'b');
