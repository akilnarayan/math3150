% Simulates the solution to the one-dimensional constant-coefficient,
% homogeneous wave equation.

clear
close all

N = 128; % Number of eigenvalue/eigenfunction pairs to compute
L = 1; % Length of the domain
csp = 2.0; % u_tt = csp^2 u_xx

% Initial data -- should be a function defined in a vectorized way
% u(x,0) = f(x)
% u_t(x,0) = g(x)
%f = @(x) -3*x.*(1-x).^2.*exp(cos(4*x));
f = @(x) x>=0.5;
g = @(x) 10*x.*(1-x);

% Boundary conditions
%  a * u(0) + b * u'(0) = v(1)
%  c * u(L) + d * u'(L) = v(2)
a = 1; b = 0;
c = 0; d = 1;
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
  gg = @(x) g(x) - c(1);
end

% Solve separation of variables eigenvalue problem
[lambda0, v0] = zero_eigenvalues(A);
if isempty(lambda0)
  zero_ev = false;
else
  zero_ev = true;
end
[lambda, v] = positive_eigenvalues(A, N-length(lambda0));
lambda = [lambda0; lambda]/L^2;
v = [v0; v]/sqrt(L);

M = 2*N; % Number of quadrature points for approximation of integrals
[x,w] = gauss_quadrature(M, 0, L);

% Project ff onto span of eigenfunctions
% ff = sum_{n=1}^N c_n phi_n(x)
ffx = ff(x);
ggx = gg(x);
ac = zeros([N 1]);
bc = zeros([N 1]);
phinx = zeros([numel(x) N]);
for n = 1:N
  if (zero_ev && (n==1))
    phinx(:,n) = v(n,1) + v(n,2) * x;
  else
    phinx(:,n) = v(n,1) * cos(sqrt(lambda(n))*x) + v(n,2) * sin(sqrt(lambda(n))*x);
  end
  ac(n) = sum(w.*ffx.*phinx(:,n));
  bc(n) = sum(w.*ggx.*phinx(:,n))/csp/sqrt(lambda(n));
end

%%%%%% Visualization options
Nviz = 7;
lineprops = {'linewidth', 3};
labelprops = {'fontsize', 16, 'fontweight', 'b', 'interpreter', 'latex'};
axesprops = {'fontsize', 16, 'fontweight', 'b'};

ymin = -3.00;
ymax = 3.00;

dt = 0.005;
T = 10;
%%%%%%

figure;
subplot(2,3,1);
set(plot(x, phinx*ac, 'r'), lineprops{:});
set(xlabel('$\mathbf{x}$'), labelprops{:});
set(ylabel('$\mathbf{u(x,0) = f(x)}$'), labelprops{:});
set(gca, axesprops{:});
set(gca, 'ylim', [ymin ymax]);

subplot(2,3,2);
set(plot(x, phinx(:,1:Nviz)*diag(ac(1:Nviz))), lineprops{:});
set(xlabel('$\mathbf{x}$'), labelprops{:});
set(ylabel('$\mathbf{a_n \phi_n(x)}$ for $\mathbf{n=1, \ldots}$'), labelprops{:});
set(gca, axesprops{:});

subplot(2,3,3);
set(plot(x, phinx(:,1:Nviz)*diag(bc(1:Nviz))), lineprops{:});
set(xlabel('$\mathbf{x}$'), labelprops{:});
set(ylabel('$\mathbf{b_n \phi_n(x)}$ for $\mathbf{n=1, \ldots}$'), labelprops{:});
set(gca, axesprops{:});

subplot(2,3,4);
usol = plot(x, phinx*(ac.*cos(csp*sqrt(lambda)*0) + bc.*sin(csp*sqrt(lambda)*0)) , 'r'); set(usol, lineprops{:});
set(xlabel('$\mathbf{x}$'), labelprops{:});
set(ylabel('$\mathbf{u(x,t)}$'), labelprops{:});
utitle = title('$\mathbf{t=0}$');
set(utitle, labelprops{:});
set(gca, axesprops{:});
set(gca, 'ylim', [ymin ymax]);

subplot(2,3,5);
cosphins = plot(x, phinx(:,1:Nviz)*diag(ac(1:Nviz)));
for j = 1:numel(cosphins); set(cosphins(j), lineprops{:}); end;
set(xlabel('$\mathbf{x}$'), labelprops{:});
set(ylabel('$\mathbf{a_n \cos(c \sqrt{\lambda_n} t) \phi_n(x)}$ for $\mathbf{n=1, \ldots}$'), labelprops{:});
set(gca, axesprops{:});

subplot(2,3,6);
sinphins = plot(x, phinx(:,1:Nviz)*diag(sin(0)*bc(1:Nviz)));
for j = 1:numel(sinphins); set(sinphins(j), lineprops{:}); end;
set(xlabel('$\mathbf{x}$'), labelprops{:});
set(ylabel('$\mathbf{b_n \sin(c \sqrt{\lambda_n} t) \phi_n(x)}$ for $\mathbf{n=1, \ldots}$'), labelprops{:});
set(gca, axesprops{:});


pause

t = 0;
while (T - t) > 0

  t = min(T, t + dt);
  subplot(2,3,4);
  ax = axis();
  set(usol, 'ydata', phinx*(ac.*cos(csp*sqrt(lambda)*t) + bc.*sin(csp*sqrt(lambda)*t)));
  set(utitle, 'string', ['$\mathbf{t=' sprintf('%1.2f', t) '}$']);
  axis(ax);

  subplot(2,3,5);
  ax = axis();
  for j = 1:numel(cosphins)
    set(cosphins(j), 'ydata', phinx(:,j)*ac(j)*cos(csp*sqrt(lambda(j))*t));
  end
  axis(ax);

  subplot(2,3,6);
  ax = axis();
  for j = 1:numel(sinphins)
    set(sinphins(j), 'ydata', phinx(:,j)*bc(j)*sin(csp*sqrt(lambda(j))*t));
  end
  axis(ax);
  drawnow;


end
