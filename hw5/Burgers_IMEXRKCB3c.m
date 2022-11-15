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
y_march = -sin(pi * x / L) - sin(2 * pi * x / L) + sin(6 * pi * x / L);
NR_PlotXY(x,y_march,0,0,L,-3,3)
% Precalculate the time-stepping coefficients used in the simulation
% h_bar = dt*[8/15 2/15  1/3];     d=h_bar/(2*dx^2);          a= -h_bar/(2*dx^2);
% beta_bar = [1    25/8  9/4];     e=beta_bar.*h_bar/(2*dx);  b=1+h_bar/dx^2;
% zeta_bar = [0   -17/8 -5/4];     f=zeta_bar.*h_bar/(2*dx);  c= -h_bar/(2*dx^2);
% if last two characters end in bt, it means it's from Butcher tableau
aexbt = [0, 0, 0, 0; ...
    3375509829940 / 42525919076317, 0, 0, 0 ; ...
    0, 272778623835 / 1039454778728, 0, 0; ...
    0, 0, 0, 0];
aimbt = [0, 0, 0, 0; ...
    0, 3375509829940 / 4525919076317, 0, 0;
    0, 11712383888607531889907 / 32694570495602105556248, 566138307881 / 912153721139, 0; ...
    0, 0, 1660544566939 / 2334033219546, 0];
bbt = [0, 673488652607 / 2334033219546, 493801219040 / 853653026979, 184814777513 / 1389668723319];
cbt = [0, 3375509829940 / 42525919076317, 272778623835 / 1039454778728, 1];
% bimbt2o = [0, 366319659506 / 1093160237145, 270096253287 / 480244073137, 104228367309 / 1017021570740];
% bexbt2o = [449556814708 / 1155810555193, 0, 210901428686 / 1400818478499, 480175564215 / 1042748212601];

exfun_denom = (dx)^2;
imfun_denom = 2 * dx;
t1 = x(2:N);
for tStep= 1:Tmax / dt
    for k = 1:4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL 4 RK SUBSTEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t0 = x(2:N);
        t2 = -y_march(2:N) .* (y_march(3:N + 1) - y_march(1:N - 1)) .* aimbt(k, k);  % nonlinear for EX
        if k > 2
            y_march(1:N - 1) = t0;
        end
    end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF RK LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_march(1:N - 1) = t1;
    if (mod(tStep,PlotInterval)==0) NR_PlotXY(x,y_march,tStep*dt,0,L,-3,3); end
end
end % function NR_Burgers_CNRKW3_FD
