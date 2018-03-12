function[lambda] = rootfinding(f, intervals, tol)
% rootfinding -- Computes roots of a transcendental equation
%
% lambda = rootfinding(f, intervals, tol)
%
%   The input f is a function handle evaluting a scalar function of one
%   variable. This function uses matlab's fzero to find N roots of this
%   equations, where N == size(intervals, 1). Here, intervals is a N x 2 matrix, where each row is a bounding interval for the root.
%
%   The last input tol is a scalar that denotes how far from the interval edges
%   to begin searching. (The assumption is that endpoints of the interval are
%   poles for f.)

N = size(intervals, 1);
assert(size(intervals, 2) == 2);

lambda = zeros([N 1]);

for n = 1:N
  lambda(n) = fzero(f, [intervals(n,1) + tol, intervals(n,2) - tol]);
end
