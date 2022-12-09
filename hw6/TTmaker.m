function A = TTmaker(a, b, c, n)
% Creates tridiagonal Toeplitz matrix 
% a for subdiagonal, b for main diagonal, c for superdiagonal
A = a * diag(ones(n - 1, 1), -1) + b * diag(ones(n, 1), 0) + c * diag(ones(n - 1, 1), 1);
end

