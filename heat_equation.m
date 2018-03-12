% Simulates the solution to the one-dimensional constant-coefficient,
% homogeneous heat equation.

clear
close all

N = 128; % Number of eigenvalue/eigenfunction pairs to compute
L = 1; % Length of the domain
k = 0.01; % u_t = k u_xx

% Initial data -- should be a function defined in a vectorized way
f = @(x) x.*(1-x).*exp(cos(4*x));

% Boundary conditions
%  a * u(0) + b * u'(0) = v(1)
%  c * u(L) + d * u'(L) = v(2)
a = 1; b = 0;
c = 1; d = 0;
v = [0; 0];

tol = 1e-12;

% Preprocessing -- normalize boundary data
temp = norm([a b]); assert (temp > tol, 'Must specify nontrivial boundary conditions');
a = a/temp; b = b/temp; v(1) = v(1)/temp;
temp = norm([c d]); assert (temp > tol, 'Must specify nontrivial boundary conditions');
c = c/temp; d = d/temp; v(2) = v(2)/temp;

A = [a b/L; c d/L];
Abc = [b a; c*L+d c];

% Compute equilibrium solution to handle non-homogenous boundary conditions
c = pinv(Abc)*v;
if norm(Abc*c - v) > 1e-12
  error('Boundary conditions do not satisfy solvability condition');
else
  ff = @(x) f(x) - c(1)*x - c(2);
end

% Solve separation of variables eigenvalue problem
[lambda0, v0] = zero_eigenvalues(A);
[lambda, v] = positive_eigenvalues(A, N-length(lambda0));
lambda = [lambda0; lambda]/L^2;
v = [v0; v]/sqrt(L);

M = 2*N; % Number of quadrature points for approximation of integrals
[x,w] = gauss_quadrature(M, 0, L);

% Project ff onto span of eigenfunctions
% ff = sum_{n=1}^N c_n phi_n(x)
ffx = ff(x);
c = zeros([N 1]);
phinx = zeros([numel(x) N]);
for n = 1:N
  phinx(:,n) = v(n,1) * cos(sqrt(lambda(n))*x) + v(n,2) * sin(sqrt(lambda(n))*x);
  c(n) = sum(w.*ffx.*phinx(:,n));
end

%%%%%% Visualization options
Nviz = 8;
lineprops = {'linewidth', 3};
labelprops = {'fontsize', 16, 'fontweight', 'b', 'interpreter', 'latex'};
axesprops = {'fontsize', 16, 'fontweight', 'b'};

dt = 0.005;
T = 10;
%%%%%%

figure;
subplot(2,2,1);
set(plot(x, phinx*c, 'r'), lineprops{:});
set(xlabel('$\mathbf{x}$'), labelprops{:});
set(ylabel('$\mathbf{u(x,0) = f(x)}$'), labelprops{:});
set(gca, axesprops{:});

subplot(2,2,2);
set(plot(x, phinx(:,1:Nviz)*diag(c(1:Nviz))), lineprops{:});
set(xlabel('$\mathbf{x}$'), labelprops{:});
set(ylabel('$\mathbf{c_n \phi_n(x)}$ for $\mathbf{n=1, \ldots}$'), labelprops{:});
set(gca, axesprops{:});

subplot(2,2,3);
usol = plot(x, phinx*c, 'r'); set(usol, lineprops{:});
set(xlabel('$\mathbf{x}$'), labelprops{:});
set(ylabel('$\mathbf{u(x,t)}$'), labelprops{:});
utitle = title('$\mathbf{t=0}$');
set(utitle, labelprops{:});
set(gca, axesprops{:});

subplot(2,2,4);
phins = plot(x, phinx(:,1:Nviz)*diag(c(1:Nviz)));
for j = 1:numel(phins); set(phins(j), lineprops{:}); end;
set(xlabel('$\mathbf{x}$'), labelprops{:});
set(ylabel('$\mathbf{c_n \exp\left(-\lambda_n^2 t\right) \phi_n(x)}$ for $\mathbf{n=1, \ldots}$'), labelprops{:});
set(gca, axesprops{:});

pause

t = 0;
while (T - t) > 0

  t = min(T, t + dt);
  subplot(2,2,3);
  ax = axis();
  set(usol, 'ydata', phinx*(c.*exp(-k*lambda*t)));
  set(utitle, 'string', ['$\mathbf{t=' sprintf('%1.2f', t) '}$']);
  axis(ax);
  drawnow;

  subplot(2,2,4);
  ax = axis();
  for j = 1:numel(phins)
    set(phins(j), 'ydata', phinx(:,j)*c(j)*exp(-k*lambda(j)*t));
  end
  axis(ax);
  drawnow;

end
