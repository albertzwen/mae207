function f = lorenz(x)
%Simulates based on the Lorenz equation 
%   Algorithm 10.4
% if varargin == 3
%     sigma = 4;
%     b = 1;
%     r = 48;
% end
sigma = 4;
b = 1;
r = 48;
x1 = x(1);
x2 = x(2);
x3 = x(3);
f = [sigma * (x2 - 1), ...
    -x2 - x1 * x3, ...
    -b * x3 + x1 * x2 - b * r];
end

