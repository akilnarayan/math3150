% Deno for approximation by Fourier Series
%
% Computes truncated Fourier series approximations for a given function f(x) on
% the interval [-L, L].

clear
close all

L = 1;
f = @(x) (x >= L/2);

% Truncated Fourier sums of the form
%
%          N                           N
%  f(x) = \sum   a_n cos(n \pi x/L) + \sum    b_n sin(n \pi x/L)
%          n=0                        n=1
%
% are computed for various values of N.

M = 1e3;
Nmax = 100;
Ns = [3 10 20 Nmax];

assert(numel(Ns) == 4);

% Exact Fourier series coefficients
n = (1:Nmax).';
an = -sin(n*pi/2)./(n*pi);
bn = -1./(n*pi) .* ( (-1).^n - cos(n*pi/2) );

an = [1/4; an]; % Concatenate with a_0


x = linspace(-L, L, M).';
fn = zeros([M 4]);
for r = 1:4

  fn(:,r) = an(1);
  for n = 1:Ns(r);
    fn(:,r) = fn(:,r) + an(n+1) * cos(n*pi*x/L) + bn(n) * sin(n*pi*x/L);
  end

end

fx = f(x);

figure; 
semilogy(1:n, abs(bn), 'b.', 'markersize', 20); hold on;
semilogy(0:n, abs(an), 'r.', 'markersize', 20); 
set(legend('$|b_n|$', '$|a_n|$'), 'interpreter', 'latex');
legend boxoff
set(gca, 'fontsize', 20, 'fontweight', 'b');
set(ylabel('Coefficient magnitude'), 'interpreter', 'latex');
set(xlabel('Coefficient index $n$'), 'interpreter', 'latex');

figure;
for r = 1:4
  subplot(2,2,r);
  plot(x, fx, 'k:', 'linewidth', 3); hold on;
  plot(x, fn(:,r), 'b', 'linewidth', 2);
  set(xlabel('$x$'), 'interpreter', 'latex');
  set(title(['$N = ' num2str(Ns(r)) '$']), 'interpreter', 'latex');
  if r == 1 
    set(legend('Function $f$', 'Fourier series'), 'interpreter', 'latex');
    legend boxoff
  end
  set(gca, 'fontsize', 16, 'fontweight', 'b');
end
