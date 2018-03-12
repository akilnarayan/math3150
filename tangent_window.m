function[left, center, right] = tangent_window(index)
% tangent_window -- Computes bounding interval for tangent function
%
% [left, center, right] = tangent_window(index)
%
%   Returns an interval and center [left, center, right] such that the tangent
%   function is bijective on that interval with range (-infty, infty). The
%   index==0 window is [-pi/2, pi/2].  These windows can be defined by:
%
%     [left, center, right] = [ (2*index - 1)*pi/2, index*pi, (2*index + 1)*pi/2]

left =  (2*index - 1)*pi/2;
center = index *pi;
right = (2*index + 1)*pi/2;
