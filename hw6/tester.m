%% hw6_att1 with Dirichlet BCs
% for k=1:4
%   switch k
%     case 1, N=128; dt=.01;
%     case 2, N= 32; dt=.01; 
%     case 3, N=128; dt=.02; 
%     case 4, N=128; dt=.04; 
%   end
%   disp(sprintf('Testing hw6_att1 with N=%d, dt=%0.4g, Dirichlet BCs',N,dt));
%   hw6_att1(4,4,N,dt, true); grid; if k<4, pause; else, disp(' '), end
% end

%% hw6_att1 with periodic BCs
% for k=1:4
%   switch k
%     case 1, N=128; dt=.01;
%     case 2, N= 32; dt=.01; 
%     case 3, N=128; dt=.02; 
%     case 4, N=128; dt=.04; 
%   end
%   disp(sprintf('Testing hw6_att1 with N=%d, dt=%0.4g, periodic BCs',N,dt));
%   hw6_att1(4,4,N,dt, false); grid; if k<4, pause; else, disp(' '), end
% end
%% hw6_att2 with Dirichlet BCs
% for k=1:4
%   switch k
%     case 1, N=128; dt=.01;
%     case 2, N= 32; dt=.01; 
%     case 3, N=128; dt=.02; 
%     case 4, N=128; dt=.04; 
%   end
%   disp(sprintf('Testing hw6_att2 with N=%d, dt=%0.4g, Dirichlet BCs',N,dt));
%   hw6_att2(4,4,N,dt, true); grid; if k<4, pause; else, disp(' '), end
% end
%% hw6_att2 with periodic BCs
for k=1:4
  switch k
    case 1, N=128; dt=.01;
    case 2, N= 32; dt=.01; 
    case 3, N=128; dt=.02; 
    case 4, N=128; dt=.04; 
  end
  disp(sprintf('Testing hw6_att2 with N=%d, dt=%0.4g, periodic BCs',N,dt));
  hw6_att2(4,4,N,dt, false); grid; if k<4, pause; else, disp(' '), end
end