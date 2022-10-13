function f = lorenz(x1, x2, x3, pStruct)
%Simulates based on the Lorenz equation 
%   Algorithm 10.4
f = [pStruct.sigma * (x2 - 1); ...
    -x2 - x1 * x3; ...
    -pstruct.b * x3 + x1 * x2 - pStruct.b * pStruct.r];
end

