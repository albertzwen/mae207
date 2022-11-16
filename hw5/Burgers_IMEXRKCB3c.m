function Burgers_IMEXRKCB3c
% function <a href="matlab:Burgers_IMEXRKCB3c">Burgers_IMEXRKCB3c</a>
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
y_march = -sin(pi * x / L) - sin(2 * pi * x / L) + sin(6 * pi * x / L);
NR_PlotXY(x,y_march,0,0,L,-3,3)
% crank up viscosity on second derivative coeff
% change time discretization
% Precalculate the time-stepping coefficients used in the simulation
% if last two characters end in bt, it means it's from Butcher tableau
bbt = [0, 673488652607 / 2334033219546, 493801219040 / 853653026979, 184814777513 / 1389668723319];
aexbt = [0, 0, 0, 0; ...
    3375509829940 / 42525919076317, 0, 0, 0 ; ...
    0, 272778623835 / 1039454778728, 0, 0; ...
    bbt];
aimbt = [0, 0, 0, 0; ...
    0, 3375509829940 / 4525919076317, 0, 0;
    0, 11712383888607531889907 / 32694570495602105556248, 566138307881 / 912153721139, 0; ...
    bbt(1), bbt(2), 1660544566939 / 2334033219546, 0];
cbt = [0, 3375509829940 / 42525919076317, 272778623835 / 1039454778728, 1];
% bimbt2o = [0, 366319659506 / 1093160237145, 270096253287 / 480244073137, 104228367309 / 1017021570740];
% bexbt2o = [449556814708 / 1155810555193, 0, 210901428686 / 1400818478499, 480175564215 / 1042748212601];
% constants from CN/RKW3
hbar = dt .* [cbt(2), cbt(3) - cbt(2), 1 - cbt(3)];
betabar = [aimbt(2, 1) / cbt(2), aimbt(3, 2) / (cbt(3) - cbt(2)), bbt(3) / (1 - cbt(3))];
zetabar = [0, -17 / 8 / (cbt(3) - cbt(2)), -5 / 4 / (1 - cbt(3))];
dxsquared = (dx)^2;
dxmult2 = 2 * dx;
atdiag = -hbar ./ (2 * dxsquared);
btdiag = 1 + hbar ./ dxsquared;
ctdiag = -hbar ./ (2 * dxsquared);

for tStep= 1:Tmax / dt
    for k = 1:4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL 4 RK SUBSTEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind = k;
        r = -y_march(2:N) .* (y_march(3:N + 1) - y_march(1:N - 1)); % nonlinear
        % QUESTIONABLE
%         atdiag = 0;
%         btdiag = 1.07;
%         ctdiag = 0;

        rhs = y_march(2:N);
        
        if (k > 2)
            ex_weight = (aimbt(k, k - 1) - bbt(k - 1)) .* dt;
            im_weight = (aexbt(k, k - 1) - bbt(k - 1)) .* dt;

%             atdiag = -ex_weight  / dxmult2;
%             btdiag = 1 + im_weight / dxsquared;
%             ctdiag = -ex_weight / dxmult2; % questionable

            rhs = rhs + ...
                im_weight .* (y_march(3:N + 1) - 2 * y_march(2:N) + y_march(1:N - 1)) + ...
                ex_weight .*   r;
        end
        if k == 4
            ind = 3;    % not sure why 
        end
        y_march(2:N) = NR_ThomasTT(atdiag(ind), btdiag(ind), ctdiag(ind), rhs', N - 1);
    end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF RK LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     y_march(2:N) = rhs;
    % enforce Dirichlet homogenous BCs
    y_march(1) = 0;
    y_march(N + 1) = 0;
    if (mod(tStep,PlotInterval)==0) NR_PlotXY(x,y_march,tStep*dt,0,L,-3,3); end
end
end % function NR_Burgers_CNRKW3_FD
