function hw6_att2(L, Tmax, N, dt, BCcase_a_bool)
% Solves the 1D wave equation using CN in time, compact Pade in space
% Now with interleaving q and v
% L = dimension width
% Tmax = total simulation time
% N = number of gridpoint discretizations
% dt = time step size
dx = L/N; 
t = 0;
xgrid = (-N / 2:N / 2 - 1)' * dx;   % discretized points on x
q=exp(-xgrid.^2/0.1); 
v = 0;
x = zeros(2 * N, 1);
% interleaving loop
diff = 0;
for i = 1:2 * N
    if mod(i, 2) == 1
        if i == 1
            x(i) = q(i);
        else
            x(i) = q(i - diff);
        end
        diff = diff + 1;
    else
        x(i) = v(i - diff);
    end
    
end
% time marching loop
for step = 1:Tmax / dt
    x = A \ r;
    t = t + dt;
    if BCcase_a_bool
        % enforce homogenous Dirichlet BCs
        x([1 N], 1) = 0;
    else
        % enforce periodic BCs
        
    end
    NR_PlotXY(xgrid,x(1:N),t,-L/2,L/2,-0.2,1.2); 
end
end

