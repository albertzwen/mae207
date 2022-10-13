function f = rossler(x1, x2, x3, a, b, c)
%Simulates based on the Rossler equation 
%   Algorithm 10.4
if vargin == 3
    a = 0.2;
    b = 0.2;
    c = 5.7;
end
f = [-x2 - x3, 
    x1 + a * x2,
    b + x3 * (x1 - c)];
end