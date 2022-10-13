function f = lorenz(x1, x2, x3, sigma, b, r)
%Simulates based on the Lorenz equation 
%   Algorithm 10.4
f = [sigma * (x2 - 1), ...
    -x2 - x1 * x3, ...
    -b * x3 + x1 * x2 - b * r];
end

