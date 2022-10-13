function f = rossler(x1, x2, x3, pStruct)
%Simulates based on the Rossler equation 
%   Algorithm 10.4
f = [-x2 - x3; x1 + pStruct.a * x2; pStruct.b + x3(x1 - pStruct.c)];
end