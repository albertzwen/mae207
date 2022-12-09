function hw6_att1(L, Tmax, N, dt, BCcase_a_bool)
% Solves the 1D wave equation using CN in time, compact Pade in space
% L = dimension width
% Tmax = total simulation time
% N = number of gridpoint discretizations
% dt = time step size
dx = L/N; 
t = 0;
xgrid = (-N / 2:N / 2 - 1)' * dx;   % discretized points on x
% Initial conditions
q=exp(-xgrid.^2/0.1); 
% v=0;    
v = zeros(N, 1);    % is this a valid IC?
x = [q; v]; % vector x, q = x(1:N)
NR_PlotXY(xgrid,x(1:N),t,-L/2,L/2,-0.2,1.2); 
C = TTmaker(0.1, 1, 0.1, N);
D = 1 / (dx)^2 * TTmaker(1.2, -2.4, 1.2, N);
A = 1 / dt * [eye(N), zeros(N); zeros(N), C] ...
    - 0.5 * [zeros(N), eye(N); D, zeros(N)];
rr = (1 / dt * [eye(N), zeros(N); zeros(N), C] ...
    + 0.5 * [zeros(N), eye(N); D, zeros(N)]);
for step = 1:Tmax / dt
    % implicit solve at each time step
    x = A \ (rr * x);   % r vector is rr * x
    t = t + dt;
    if BCcase_a_bool
        % enforce homogenous Dirichlet BCs
        x([1 N], 1) = 0;    
%     else
        % enforce periodic BCs- How?  extended diagonal, circulant
        
    end
    NR_PlotXY(xgrid,x(1:N),t,-L/2,L/2,-0.2,1.2); 
end
end

