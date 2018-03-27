function[lambda, k] = positive_eigenvalues(A, N)
% positive_eigenvalues -- Determines positive eigenvalues of an ODE
%
% [lambda, k] = positive_eigenvalues(A, N)
%
%   Determines the first N eigenvalues lambda > 0 that are eigenvalues for the
%   boundary value problem
%
%          y''(x) + lambda*y(x) = 0,   0 < x < 1
%
%           A(1,1)*y(0) + A(1,2)*y'(0) = 0,
%           A(2,1)*y(1) + A(2,2)*y'(1) = 0.
%
%   The output k is an N x 2 matrix, such that 
%
%     y(x) = k(n,1) cos( sqrt(lambda(n)) * x) + k(n,2) * sin( sqrt(lambda(n)) * x)
%
%   is an orthonormal eigenfunction corresponding to eigenvalue lambda(n) for n = 1, ..., N

lambda = zeros([N 1]);
k = zeros([N 2]);

[a,b,c,d] = deal( A(1,1), A(1,2), A(2,1), A(2,2) );

% Temporary variables used by rootfinders
tol = 1e-13;
intervals = zeros([N 2]);

if ( (a == 0) && (c == 0) ) || ( (b == 0) && (d == 0) )

  lambda = ( (1:N).' * pi ).^2;

elseif ( (a == 0) && (d == 0) ) || ( (b == 0) && (c == 0) )

  lambda = ( (2*(1:N).'-1) * pi/2 ).^2;

elseif ( (a == 0) || (c == 0) )

  % Solve x*tan(x) = (c/d - a/b)
  f = @(x) x*tan(x) - (c/d - a/b);

  % There is an n=1 solution for 0 < x < pi/2 iff (c/d - a/b) > 0
  if (c/d - a/b) > 0
    [left, center, right] = tangent_window((0:N-1)');
    intervals = [center, right];
  else
    [left, center, right] = tangent_window((1:N)');
    intervals = [left center];
  end

  lambda = rootfinding(f, intervals, tol).^2;

elseif ( (b == 0) || (d == 0) )

  % Solve tan(x)/x = (b/a - d/c)
  f = @(x) tan(x)./x - (b/a - d/c);

  % There is an n=1 solution for 0 < x < pi/2 iff (b/a - d/c) > 1
  if (b/a - d/c) > 1
    [left, center, right] = tangent_window((0:(N-1))');
    intervals(1,:) = [center(1), right(1)];
    intervals(2:N,:) = [left(2:N), right(2:N)];
  else
    [left, center, right] = tangent_window((1:N)');
    intervals = [left right];
  end

  lambda = rootfinding(f, intervals, tol).^2;

else % Then a, b, c, d all do not vanish
  % Solve tan(x)/x = (bc - ad) / (ac + bd*x^2)
  f = @(x) tan(x)./x - (b*c - a*d)./(a*c + b*d*x.^2);

  if (a*c)/(b*d) > 0 % Then whole thing is always single-signed

    % There is an n=1 solution for 0 < x < pi/2 if this condition is satisfied
    %if (b*c - a*d)./(a*c + b*d) > 0
    if (b*c - a*d)./(a*c) > 1
      [left, center, right] = tangent_window((0:(N-1))');
      intervals(1,:) = [center(1) right(1)];
      intervals(2:N,:) = [left(2:N) right(2:N)];
    else
      [left, center, right] = tangent_window((1:N)');
      intervals = [left right];
    end

    lambda = rootfinding(f, intervals, tol).^2;

  else
    % Pole at x = sqrt( -ac/bd)
    pole = sqrt(-a*c/(b*d));

    % Which tangent window is this pole inside?
    index = ceil((pole*2/pi-1)/2);
    upper = pi/2 + pi*index;
    if index > 0
      lower = upper - pi;
    else
      lower = 0;
    end

    % At x=0, RHS is (b/a-d/c)
    % At x=infty, RHS is 1/x^2 * (c/d-a/b)
    % When -ac/(bd) > 0, then 
    % (b/a-d/c) < 0 =====> (c/d-a/b) > 0, and vice versa

    if (b/a - d/c) < 0 
      % - no pole in window 0
      % - 1 pole in windows 1, ..., index-1
      % - 2 poles in window index
      % - 1 pole in window index+1, ...

      M = index-1;
      [left, center, right] = tangent_window((1:M)');
      intervals(1:M,:) = [left right];

      intervals(M+1,:) = [left(end) + pi pole];
      intervals(M+2,:) = [pole upper];

      [left, center, right] = tangent_window(((M+2):(N-1))');
      intervals((M+3):end, :) = [left right];

      intervals = intervals(1:N,:);

    else % (b/a - d/c) > 0 
      ev_count = 0;

      if index == 0
        % Suppose index == 0
        % - Test f(0+tol) vs f(pole-tol) for solution
        % - 1 pole in window index+1, ...

        if (b/a - d/c) < 1 % Then an intersection happens in window 0
          intervals(1,:) = [0 pole];
          [left, center, right] = tangent_window((1:(N-1))');
          intervals(2:N,:) = [left right];
        else
          [left, center, right] = tangent_window((1:N)');
          intervals(1:N,:) = [left right];
        end
          
      elseif (b/a - d/c) > 1 % And also index > 0
        % Suppose index > 0, (b/a- d/c) > 1
        % - 1 pole in window 0, 1, ..., index-1
        % - 0 pole in window index
        % - 1 pole in window index+1, ...

        if index < N
          [left, center, right] = tangent_window([(0:(index-1))'; ((index+1):N)']);
        else
          [left, center, right] = tangent_window([(0:(index-1))']);
        end
        intervals(1,:) = [center(1), right(1)];
        intervals(2:N,:) = [left(2:N) right(2:N)];

      else
        % Suppose index > 0, (b/a - d/c) < 1
        % - 0 pole in window 0
        % - 1 pole in window 1, ..., index-1
        % - 0 pole in window index
        % - 1 pole in window index+1, ...

        [left, center, right] = tangent_window([(1:(index-1))'; ((index+1):(N+1))']);
        intervals(1:N,:) = [left(1:N) right(1:N)];

      end

    end

    lambda = rootfinding(f, intervals, tol).^2;

  end
end

% Compute eigenfunctions
nu = sqrt(lambda);
k = [ b*nu -a*ones([N 1]) ];
knorm = 1/2 * (b^2*lambda+ a^2 + sinc(2/pi * nu).*(b^2.*lambda - a^2) + a*b*(cos(2*nu) - 1));

k(:,2) = k(:,2)./sqrt(knorm);
k(:,1) = k(:,1)./sqrt(knorm);
