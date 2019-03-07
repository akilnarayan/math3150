% Video demonstration of the heat equation solution on an infinite domain.
%
% The solution is computed as the convolution of the initial data against the
% heat kernel.
%
% The equation is
%
%   u_t = k * u_xx,             u(x,0) = f(x)
%
% for -\infty < x < \infty, and t > 0. k is a positive constant
%
% Uses a discretization method to solve over x in [-2,2].

clear
close all

k = 1;
%f = @(x) cos(6*pi*x) .* (abs(x) <= 1/12);
%f = @(x) (1 + cos(6*pi*x)) .* (abs(x) <= 6/12);
f = @(x) abs(x) <= 0.2;


h = @(x,t) 1/sqrt(4*k*pi*t) * exp(-x.^2./(4*k*t));

[s,w] = hermite_gauss_quadrature(500);

x = linspace(-2, 2, 301).';
u = double(f(x));

T = 1;
dt = 0.001;
t = 0;

subplot(3, 1, 1);
plot(x, u, 'r', 'linewidth', 2); 
set(xlabel('$x$'), 'interpreter', 'latex');
set(title('$u(x,0) = f(x)$'), 'interpreter', 'latex');
set(gca, 'fontsize', 16, 'fontweight', 'b');

subplot(3, 1, 2);
h0 = zeros(size(x)); h0(x==0) = 4;
hplot = plot(x, h0, 'b', 'linewidth', 2); 
set(xlabel('$x$'), 'interpreter', 'latex');
time = sprintf('%1.4f', t);
set(title(['$h(x,t=' time ')$']), 'interpreter', 'latex');
set(gca, 'fontsize', 16, 'fontweight', 'b');
temp = axis;
hymin = temp(3);
hymax = temp(4);

subplot(3, 1, 3);
uplot = plot(x, u, 'k', 'linewidth', 2);
set(xlabel('$x$'), 'interpreter', 'latex');
set(gca, 'fontsize', 16, 'fontweight', 'b');
time = sprintf('%1.4f', t);
set(title(['$u(x,t=' time ')$']), 'interpreter', 'latex');
temp = axis;
uymin = temp(3);
uymax = temp(4);
pause

while t < T
  t = t + dt;

  for qx = 1:length(x)
    u(qx) = sum(w.*f(x(qx) - 2*s*sqrt(k*t)))/sqrt(pi);
  end

  set(uplot, 'ydata', u);
  set(hplot, 'ydata', h(x,t));

  subplot(3,1,2);
  temp = axis;
  axis([temp(1:2) hymin hymax]);
  time = sprintf('%1.4f', t);
  set(title(['$h(x,t=' time ')$']), 'interpreter', 'latex');

  subplot(3,1,3);
  temp = axis;
  axis([temp(1:2) uymin uymax]);
  time = sprintf('%1.4f', t);
  set(title(['$u(x,t=' time ')$']), 'interpreter', 'latex');
  drawnow;

end
