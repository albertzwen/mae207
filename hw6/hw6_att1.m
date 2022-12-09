function Wave1D_CN_Pade(L, Tmax, c, N, dt)
% L = dimension width
% Tmax = total simulation time
% N = number of gridpoint discretizations
% dt = time step size
dx = L/N; 
IterSteps = 2;
t = 0;
x = (-N/2:N/2-1)'*dx; 
q=exp(-x.^2/0.1); v=0;
NR_PlotXY(x,q,t,-L/2,L/2,-0.2,1.2); 
vs=v; 
qs=v; 
a=0.6*dt*c^2/dx^2; 
b=-1.2*dt*c^2/dx^2;
for 
end

