function Burgers_IMEXRKCB3c
% Simulate the 1D Burgers on 0<x<L with homogeneous Dirichlet BCs using IMEXRKCB3c in time
% (explicit on nonlinear terms, implicit on linear terms)

% Initialize the simulation parameters (user input)
L = 100;
Tmax = 50;
N = 100;
dt = 0.5;
PlotInterval=1;
dx = L / N;
x = (0:N) .* dx; % length N + 1
% STEP 1: Discretization of unknown variable on spatial grid
u = -sin(pi * x / L) - sin(2 * pi * x / L) + sin(6 * pi * x / L);
figure;
NR_PlotXY(x,u,0,0,L,-3,3)
% TIPS
% crank up viscosity on second derivative coeff
% change time discretization

% Precalculate the time-stepping coefficients used in the simulation
% Butcher tableau of IMEXRKCB3c from CB15.pdf
% Last two characters bt => Butcher tableau
bbt= [0, 673488652607 / 2334033219546, 493801219040 / 853653026979, 184814777513 / 1389668723319];
cbt = [0, 3375509829940 / 42525919076317, 272778623835 / 1039454778728, 1];
a_exbt = [0, 0, 0, 0; ...
    cbt(2), 0, 0, 0 ; ...
    0, cbt(3), 0, 0; ...
    bbt];
a_imbt = [0, 0, 0, 0; ...
    0, 3375509829940 / 4525919076317, 0, 0;
    0, 11712383888607531889907 / 32694570495602105556248, 566138307881 / 912153721139, 0; ...
    bbt(1), bbt(2), 1660544566939 / 2334033219546, 0];
% derivations of constants through CN of RK4
zeta = [0, bbt(1) - a_exbt(2, 1), bbt(2) - a_exbt(3, 2), bbt(3) - a_exbt(4, 3)];
h_bar = dt .* [cbt(2), cbt(3) - cbt(2), cbt(4) - cbt(3), 1 - cbt(4)];
beta_bar = [a_exbt(2, 1) / cbt(2), a_exbt(3, 2) / (cbt(3) - cbt(2)), ...
    a_exbt(4, 3) / (cbt(4) - cbt(3)), bbt(4) / (1 - cbt(4))];
zeta_bar = [0, zeta(2) / (cbt(3) - cbt(2)), zeta(3) / (cbt(4) - cbt(3)), ...
    zeta(4) / (1 - cbt(4))];

dxsquared = (dx)^2;
dxmult2 = 2 * dx;

% d = h_bar ./ (2 * dxsquared);
% e = beta_bar .* h_bar ./ (dxmult2);
% f = zeta_bar .* h_bar ./ (dxmult2);
% a = -h_bar ./ (2 * dxsquared);
% b = 1 + h_bar ./ dxsquared;
% c = -h_bar / (2 * dxsquared);

y = zeros(size(x));
z = zeros(size(x));
% procedures for implementation from equation 19 of CB15.pdf
for tStep= 1:Tmax / dt
    for k = 1:4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL 4 RK SUBSTEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        r = -u(2:N) .* (u(3:N + 1) - u(1:N - 1));   % nonlinear term
        if k == 1 % register 1
%             y(2:N) = u(2:N) + d(k) .* (u(3:N + 1) - 2 .* u(2:N) + u(1:N - 1)) + e(k) .* r;
            y(2:N) = u(2:N);
        else    % register 2
            y(2:N) = u(2:N) + (a_imbt(k, k - 1) - bbt(k - 1)) .* dt .* ...
                (u(3:N + 1) - 2 .* u(2:N) + u(1:N - 1)) + ...
                (a_exbt(k, k - 1) - bbt(k - 1)) .* dt .* r;
            %             y(2:N) = u(2:N) + d(k) .* ...
            %                 (u(3:N + 1) - 2 .* u(2:N) + u(1:N - 1)) + ...
            %                 e(k).* r + f(k) .* y(2:N);
        end

        %         u(2:N) = (NR_ThomasTT(a(k), b(k), c(k), y(2:N)', N - 1))';
        z(2:N) = (NR_ThomasTT(a_imbt(k,k) * dt / (2 * dxsquared), ...
            1 - a_imbt(k,k) * dt / dxsquared, a_imbt(k, k) * dt / (2 * dxsquared), ...
            y(2:N)', N - 1))' .* ...
            y(2:N);
        y(2:N) = -(y(2:N) + a_imbt(k, k) .* z(2:N)) .* ...
            (y(3:N + 1) + a_imbt(k, k) .* (z(3:N + 1)) - ...
            (y(1:N - 1) + a_imbt(k, k) .* z(1:N - 1)));
%         y(2:N) = -(y(2:N) .* (y(3:N + 1) - y(1:N - 1)));
        
        u(2:N) = u(2:N) + bbt(k) * dt .* z(2:N) + bbt(k) * dt .* y(2:N);
    end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF RK LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % enforce Dirichlet homogenous BCs
    %     u(1) = 0;
    %     u(N + 1) = 0;
    % plot
    if (mod(tStep,PlotInterval)==0) 
%         figure(tStep);
        NR_PlotXY(x,u,tStep*dt,0,L,-3,3); 
    end
end
end
