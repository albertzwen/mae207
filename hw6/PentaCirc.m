function [A] = PentaCirc(a,b,c,d,e,G,n)
% function [G,a,b,c,d,e] = Penta(a,b,c,d,e,G,n)
% Constructs a pentadiagonal Toeplitz circulant matrix
A=diag(a(3:n), -2) + diag(b(2:n),-1)+diag(c,0)+diag(d(1:n-1),1) + diag(e(1:n - 2), 2); 
A(1,n)=a(1); A(n,1)=e(n);
end 
