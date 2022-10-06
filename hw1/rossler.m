function f = rossler(x1, x2, x3, pStruct)
%Simulates based on the Lorenz equation 
%   Needs work
f = [-x2 - x3; x1 + pStruct.a * x2; pStruct.b + x3(x1 - pStruct.c)];
end