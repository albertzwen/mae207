function hw6_att2(L, Tmax, N, dt, BCcase_a_bool)
% Solves the 1D wave equation using CN in time, compact Pade in space
% Now with interleaving q and v
% L = dimension width
% Tmax = total simulation time
% N = number of gridpoint discretizations
% dt = time step size
dx = L / N; 
t = 0;
xgrid = (-N / 2:N / 2 - 1)' * dx;   % discretized points on x
q=exp(-xgrid.^2/0.1); 

% interleaving x 
x = zeros(2 * N, 1);
x(1:2:end) = q;
NR_PlotXY(xgrid,x(1:2:end),t,-L/2,L/2,-0.2,1.2); 
% Calculating E and F
F = 1.2 / dx^2 * diag(ones(2 * N - 3, 1), -3) ...   % subsubsubdiag
    - 2.4 / dx^2 * diag(rem(1:2 * N - 1, 2), -1) ...   % subdiag
    + (diag(rem(1:2 * N - 1, 2),  1) ...
    + 1.2 / dx^2 * diag(mod(0:2 * N - 1 - 1, 2), 1));  % superdiag
E = 0.1 * diag(mod(0: 2 * N - 2 - 1, 2), -2) ... % subsubdiag
    + diag(ones(2 * N, 1), 0) ... % main diag
    + 0.1 * diag(mod(0: 2 * N - 2 - 1 , 2), 2);   % supersuperdiag
% if ~BCcase_a_bool
%     % perform modifications on corners of E
%     E(1, end) = ?;
%     E(end, 1) = ?;
% end
A = E ./ dt - F ./ 2;
rr = E ./ dt + F ./ 2;
% time marching loop
for step = 1:Tmax / dt
    x = A \ (rr * x);
    t = t + dt;
%     if BCcase_a_bool
        % enforce homogenous Dirichlet BCs
%         x([1 end - 1], 1) = 0;
%     end
    NR_PlotXY(xgrid,x(1:2:end),t,-L/2,L/2,-0.2,1.2); 
end
end

