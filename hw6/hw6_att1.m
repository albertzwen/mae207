function hw6_att1(L, Tmax, N, dt, BCcase_a_bool, interval)
% Solves the 1D wave equation using CN in time, compact Pade in space
% L = dimension width
% Tmax = total simulation time
% N = number of gridpoint discretizations
% dt = time step size
% BCcase_a_bool = boolean that indicates case a or b
%% Constant definition
dx = L/N; 
t = 0;
xgrid = (-N / 2:N / 2 - 1)' * dx;   % discretized points on x
% Initial conditions
q=exp(-xgrid.^2/0.1); 
% v=0;    
v = zeros(N, 1);    
x = [q; v]; % vector x, q = x(1:N)
%% figure saving
scenarioName = strcat('plots/att1/');
if BCcase_a_bool
    scenarioName = strcat(scenarioName, 'case_a');
else
    scenarioName = strcat(scenarioName, 'case_b');
end
scenarioName = strcat(scenarioName, '/N', num2str(N), 'dt', num2str(dt));
NR_PlotXY(xgrid,x(1:N),t,-L/2,L/2,-0.2,1.2); 
saveas(gcf, strcat(scenarioName, '/t', num2str(t), '.png'));

C = TTmaker(0.1, 1, 0.1, N);
D = 1 / (dx)^2 * TTmaker(1.2, -2.4, 1.2, N);
if ~BCcase_a_bool
    C(1, end) = 0.1;
    C(end, 1) = 0.1;
    D(1, end) = 1.2 / dx^2;
    D(end, 1) = 1.2 / dx^2;
end
A = 1 / dt * [eye(N), zeros(N); zeros(N), C] ...
    - 0.5 * [zeros(N), eye(N); D, zeros(N)];    % what does it look like?
rr = (1 / dt * [eye(N), zeros(N); zeros(N), C] ...
    + 0.5 * [zeros(N), eye(N); D, zeros(N)]);
for step = 1:Tmax / dt
    % implicit solve at each time step
    x = A \ (rr * x);   % r vector is rr * x
    t = t + dt;
    if BCcase_a_bool
        % enforce homogenous Dirichlet BCs
        x([1 N], 1) = 0;            
    end
    if mod(step, interval) == 0
        NR_PlotXY(xgrid,x(1:N),t,-L/2,L/2,-0.2,1.2); 
        saveas(gcf, strcat(scenarioName, '/t', num2str(t), '.png'));
    end
end
end

