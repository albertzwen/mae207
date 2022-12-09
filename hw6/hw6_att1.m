function hw6_att1(L, Tmax, N, dt, BCcase_a_bool)
% Solves the 1D wave equation using CN in time, compact Pade in space
% L = dimension width
% Tmax = total simulation time
% N = number of gridpoint discretizations
% dt = time step size
dx = L/N; 
t = 0;
xgrid = (-N / 2:N / 2 - 1)' * dx;   % discretized points on x
q=exp(-xgrid.^2/0.1); 
v=0;    % v = dqdt - blanking on where to get this again
x = [q; v]; % vector x 
NR_PlotXY(xgrid,x(1:N),t,-L/2,L/2,-0.2,1.2); 
C = TTmaker(0.1, 0.1, 0.1, 2 * N);
D = 1 / (dx)^2 * TTmaker(1.2, -2.4, 1.2, 2 * N);
A = 1 / dt * [eye(N), zeros(N); zeros(N), C] - 0.5 * [zeros(N), eye(N); D zeros(N)];
r = 1 / dt * [eye(N), zeros(N); zeros(N), C] + 0.5 * [zeros(N), eye(N); D zeros(N)];
for step = 1:Tmax / dt
    x = A \ r;
    t = t + dt;
    if BCcase_a_bool
        % enforce homogenous Dirichlet BCs
        x([1 N]) = 0;
    else
        % enforce periodic BCs
        
    end
    NR_PlotXY(xgrid,x(1:N),t,-L/2,L/2,-0.2,1.2); 
end
end

