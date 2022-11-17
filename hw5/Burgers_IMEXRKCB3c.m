function Burgers_IMEXRKCB3c
% function <a href="matlab:Burgers_IMEXRKCB3c">Burgers_IMEXRKCB3c</a>
% Simulate the 1D Burgers on 0<x<L with homogeneous Dirichlet BCs using CN/RKW3 in time

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

cbt = [0, 3375509829940 / 42525919076317, 272778623835 / 1039454778728, 1];
% bimbt2o = [0, 366319659506 / 1093160237145, 270096253287 / 480244073137, 104228367309 / 1017021570740];
% bexbt2o = [449556814708 / 1155810555193, 0, 210901428686 / 1400818478499, 480175564215 / 1042748212601];
hbar = dt .* [cbt(2), cbt(3) - cbt(2), 1 - cbt(3)];
betabar = [aimbt(2, 1) / cbt(2), aimbt(3, 2) / (cbt(3) - cbt(2)), bbt(3) / (1 - cbt(3))];
zetabar = [0, -17 / 8 / (cbt(3) - cbt(2)), -5 / 4 / (1 - cbt(3))];
dxsquared = (dx)^2;
dxmult2 = 2 * dx;
atdiag = -hbar ./ dxmult2;
btdiag = 1 + hbar ./ dxsquared;
ctdiag = -hbar ./ dxmult2;

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
        
        if (k > 2)
            ex_weight = (aimbt(k, k - 1) - bbt(k - 1)) .* dt;
            im_weight = (aexbt(k, k - 1) - bbt(k - 1)) .* dt;

%             atdiag = -ex_weight  / dxmult2;
%             btdiag = 1 + im_weight / dxsquared;
%             ctdiag = -ex_weight / dxmult2; % questionable

            rhs = rhs + ...
                ex_weight .* (y_march(3:N + 1) - 2 * y_march(2:N) + y_march(1:N - 1)) + ...
                im_weight .* r;
        end
        if k == 4
            ind = 3;    % not sure why 
        end
        y_march(2:N) = NR_ThomasTT(atdiag(ind), btdiag(ind), ctdiag(ind), rhs', N - 1);
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
