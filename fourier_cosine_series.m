% Deno for approximation by Fourier Cosine series
%
% Computes truncated Fourier cosine series approximations for a given function f(x) on
% the interval [0, L].

clear
close all

L = 1;
f = @(x) x;

% Truncated Fourier sums of the form
%
%          N                       
%  f(x) = \sum   A_n cos(n \pi x/L) 
%          n=0                    
%
% are computed for various values of N.

M = 1e3;
Nmax = 100;
Ns = [3 10 20 Nmax];

assert(numel(Ns) == 4);

% Exact Fourier series coefficients
n = (1:Nmax).';
An = 2./(n*pi).^2 .* ( (-1).^n - 1 );
An = [1/2; An]; % Concatenate with A_0

x1 = linspace(0, L, M).';       % For [0, L] plots
x2 = linspace(-L, L, M).';      % For [-L, L] plots
x3 = linspace(-4*L, 4*L, M).';  % For [-4L, 4L] plots

fn1 = zeros([M 4]);
fn2 = zeros([M 4]);
fn3 = zeros([M 4]);
for r = 1:4

  fn1(:,r) = An(1);
  fn2(:,r) = An(1);
  fn3(:,r) = An(1);
  for n = 1:Ns(r);
    fn1(:,r) = fn1(:,r) + An(n+1) * cos(n*pi*x1/L);
    fn2(:,r) = fn2(:,r) + An(n+1) * cos(n*pi*x2/L);
    fn3(:,r) = fn3(:,r) + An(n+1) * cos(n*pi*x3/L);
  end

end

fx1 = f(x1);
fx2 = f(x2);
fx3 = f(x3);

figure; 
semilogy(0:n, abs(An), 'r.', 'markersize', 20); 
set(legend('$|A_n|$'), 'interpreter', 'latex');
legend boxoff
set(gca, 'fontsize', 20, 'fontweight', 'b');
set(ylabel('Coefficient magnitude'), 'interpreter', 'latex');
set(xlabel('Coefficient index $n$'), 'interpreter', 'latex');

% Generate plot over [0, L]
figure;
for r = 1:4
  subplot(2,2,r);
  plot(x1, fx1, 'k:', 'linewidth', 3); hold on;
  plot(x1, fn1(:,r), 'b', 'linewidth', 2);
  set(xlabel('$x$'), 'interpreter', 'latex');
  set(title(['$N = ' num2str(Ns(r)) '$']), 'interpreter', 'latex');
  if r == 1 
    set(legend('Function $f$', 'Fourier cosine series on $[0, L]$'), 'interpreter', 'latex');
    legend boxoff
  end
  set(gca, 'fontsize', 16, 'fontweight', 'b');
end

% Generate plot over [-L, L]
figure;
for r = 1:4
  subplot(2,2,r);
  plot(x1, fx1, 'k:', 'linewidth', 3); hold on;
  plot(x2, fn2(:,r), 'b', 'linewidth', 2);
  set(xlabel('$x$'), 'interpreter', 'latex');
  set(title(['$N = ' num2str(Ns(r)) '$']), 'interpreter', 'latex');
  if r == 1 
    set(legend('Function $f$', 'Fourier cosine series on $[-L, L]$'), 'interpreter', 'latex');
    legend boxoff
  end
  set(gca, 'fontsize', 16, 'fontweight', 'b');
end

% Generate plot over [-4L, 4L]
figure;
for r = 1:4
  subplot(2,2,r);
  plot(x1, fx1, 'k:', 'linewidth', 3); hold on;
  plot(x3, fn3(:,r), 'b', 'linewidth', 2);
  set(xlabel('$x$'), 'interpreter', 'latex');
  set(title(['$N = ' num2str(Ns(r)) '$']), 'interpreter', 'latex');
  if r == 1 
    set(legend('Function $f$', 'Fourier cosine series on $[-4L, 4L]$'), 'interpreter', 'latex');
    legend boxoff
  end
  set(gca, 'fontsize', 16, 'fontweight', 'b');
end
