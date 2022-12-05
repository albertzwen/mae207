% Illustrates how the code PentaCirc may be used to construct pentadiagonal
% and pentadiagonal circulant matrices.
% Numerical Renaissance codebase, Chapter 1, https://github.com/tbewley/NR
% Copyright 2021 by Thomas Bewley, distributed under BSD 3-Clause License. 

n=7; a=randn(n,1); b=randn(n,1); c=randn(n,1); d = rand(n, 1); e = rand(n, 1);
A_pentadiagonal = Penta([0; a(2:n)], b, [c(1:n-1); 0]) 
A_pentadiagonal_circulant=NR_Tridiag(a,b,c)

