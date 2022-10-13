function f = rossler(x)
%Simulates based on the Rossler equation 
%   Algorithm 10.4
% if varargin == 3
%     a = 0.2;
%     b = 0.2;
%     c = 5.7;
% end
a = 0.2;
b = 0.2;
c = 5.7;
x1 = x(1);
x2 = x(2);
x3 = x(3);
f = [-x2 - x3, ...
    x1 + a * x2, ...
    b + x3 * (x1 - c)];
end