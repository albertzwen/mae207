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

%% interleaving x 
x = zeros(2 * N, 1);
x(1:2:end) = q;
NR_PlotXY(xgrid,x(1:2:end),t,-L/2,L/2,-0.2,1.2); 
%% Computing E and F
F = 1.2 / dx^2 * diag(ones(2 * N - 3, 1), -3) ...
    - 2.4 / dx^2 * diag(mod(1:2 * N - 1, 2), -1) ...
    + diag(mod(1:2 * N - 1, 2), 1) + 1.2 /dx^2 * diag(mod(0:2 * N - 1 - 1, 2), 1);
E = 0.1 * diag(mod(0:2 * N - 2 - 1, 2), -2) ...
    + diag(ones(2 * N, 1), 0) ...
    + 0.1 * diag(mod(0:2 * N - 2 - 1, 2), 2);
A = (1 / dt * E) - (0.5 * F);
rr = (1 / dt * E) + (0.5 * F);
for step = 1:Tmax / dt
    x = A \ (rr * x);
    t = t + dt;
    if BCcase_a_bool
        x([1 2*N - 1]) = 0;
    end
    NR_PlotXY(xgrid,x(1:2:end),t,-L/2,L/2,-0.2,1.2); 
end
end

