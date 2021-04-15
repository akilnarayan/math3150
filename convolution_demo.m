% Video demonstration of convolution. Given functions f and g, the convolution
% of these two function is
%
% (f*g)(x) = 1/(2pi) int_{-\infty}^{\infty} f(s) g(x-s) ds.
%
% The visualization is done by taking g as the "averaging" function, and f as
% the "target" function.
%
% This function truncates the domain of integration to [-1,1], and visualizes
% on [-2,2]. This effectively assumes that the support of each function is on
% [-1,1] and performs an exact convolution.

clear
close all

f = @(s) (1 - abs(s));

g = @(s) 1/(0.4)*(abs(s) < 0.2);

t = .1;
g = @(s) 1/sqrt(4 * pi * t) * exp(-s.^2./(4*t));

% Truncates functions outside [-1,1]
fwindow = @(s) f(s).*(abs(s) <= 1);
gwindow = @(s) g(s);

[s,w] = gauss_quadrature(501, -5, 5);

x = linspace(-2, 2, 301).';
fintegrand = fwindow(s);
gintegrand = gwindow(-2 - s);

h = zeros(size(x));

figure;
subplot(2,1,1);
fplot = plot(s, fintegrand, 'b', 'linewidth', 2);
hold on;
gplot = plot(s, gintegrand, 'r', 'linewidth', 2);
axis([-2, 2, 0, 2.7]);
set(gca, 'fontsize', 16, 'fontweight', 'b');
set(xlabel('$s$'), 'interpreter', 'latex');
set(title('Integrand'), 'interpreter', 'latex');
set(legend('$f(s)$', '$g(x-s)$'), 'interpreter', 'latex');

xticks = [-2, -1, 0, 1, 2];
xticklabels = {'-2','-1','0','1','2', 'x'};
set(gca, 'xtick', xticks);

subplot(2,1,2);
hold on;
hplot = plot(-2, 0, 'k', 'linewidth', 2);
dotplot = plot(-2, 0, 'k.', 'markersize', 20);
axis([-2, 2, 0, 1]);
set(gca, 'fontsize', 16, 'fontweight', 'b');
set(xlabel('$x$'), 'interpreter', 'latex');
set(title('Convolution $f \ast g$'), 'interpreter', 'latex');
set(ylabel('$2\pi (f \ast g)$'), 'interpreter', 'latex');
pause

s2 = [s.', flipud(s).'];
prevfill = [];

for qx = 1:numel(x)

  gintegrand = gwindow(x(qx) - s);
  h(qx) = sum(gintegrand.*fintegrand.*w);

  inbetween = [zeros(size(s)).', min(flipud([fintegrand, gintegrand]), [], 2).'];

  if any(xticks == x(qx));
    newticks = xticks;
    newticklabels = xticklabels(1:5);
  else
    newticks = [-2, -1, 0, 1, 2, x(qx)];
    [newticks, order] = sort(newticks);
    newticklabels = xticklabels(order);
  end

  subplot(2,1,1);
  set(gca, 'xtick', newticks, 'xticklabel', newticklabels);
  set(gplot, 'ydata', gintegrand);
  if not(isempty(prevfill))
    set(prevfill, 'Visible', 'off');
  end
  prevfill = fill(s2, inbetween, 'k');
  prevfill.Annotation.LegendInformation.IconDisplayStyle = 'off';
  alpha(prevfill, 0.5);

  set(hplot, 'xdata', x(1:qx), 'ydata', h(1:qx));
  set(dotplot, 'xdata', x(qx), 'ydata', h(qx));
  drawnow;
  pause(0.01);

end
