function [qs, energy] = NR_SolarSystemSimulator(method, Tmax, kmax) %%%%%%%
... % Define a bunch of constants and initial conditions here
h = Tmax / kmax; 
for k = 1:kmax
    t = k * h;
    if method == 'SI4'
        for ss = 1:4
            q = q + c(ss) * h * dqdt(p, m); 
            if ss < 4
                p = p + d(ss) * h * dpdt(q, m, G);
            end
        end 
    elseif method == 'RK4'
        k1q = dqdt(p, m);
        k1p = dpdt(q, m, G);
        k2q = dqdt(p + (h / 2) * k1p, m);
        k2p = dpdt(q + (h / 2) * k1q, m, G);
        k3q = dqdt(p + (h / 2) * k2p, m);
        k3p = dpdt(q + (h / 2) * k2q, m, G);
        k4q = dqdt(p + h * k3p, m, m); 
        k4p = dpdt(q + h * k3q, m , G);
        q = q + h * (k1q / 6 + (k2q + k3q) / 3 + k4q / 6);
        p = p + h * (k1p / 6 + (k2p + k3p) / 3 + k4p / 6);
    end
    ... % Do some plotting stuff here
end
... % Do some clean up stuff here
end % function NR_SolarSystemSimulator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = dqdt(p, m)
for i = 1:9
    for j = 1:3
        x(i, j) = p(i, j) / m(i);
    end
end
end % function dqdt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = dpdt(q, m, G)
    x = zeros(9, 3);
    for i = 1:9
        for j = 1:3
            for k = 1:9
                if k ~= i
                    x(i, j) = x(i, j) + G * m(i) * m(k) * q(k, j) - ...
                        q(i, j) / norm(q(k, :)') - q(i, :)') ^ 3;
                end
            end
        end
    end
end % function dpdt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
