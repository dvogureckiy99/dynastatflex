addpath("C:\dvogureckiy99\2_Education\ITMO\DMS\SPACAR")

nodes = [linspace(0,L,N)',zeros(N,1),zeros(N,1)];
n = size(nodes,1);
elements = [];
for i=1:n-1
    elements = [elements;i,i+1];
end
nprops(1).fix   = true; 

if exist('M','var')
    nprops(n).moment = M;
elseif exist('Fx','var') && exist('Fy','var')
    nprops(int32((l_Fext/L)*N)).force = [Fx,Fy,0];
elseif exist('Stroke','var')
    nprops(n).rot_z = Stroke;
end

%flexures
eprops(1).elems = linspace(1,n-1,n-1);
eprops(1).flex = [3 4];
% eprops(1).nbeams = 1;
eprops(1).emod = double(E);
eprops(1).smod = double(G);
eprops(1).dens = double(rho);
eprops(1).cshape ='rect';
eprops(1).dim = [double(H) double(W)];
eprops(1).orien = [0 0 1];
eprops(1).color = 'blue';
eprops(1).cw = true;

opt.filename = 'simple_beam';
opt.loadsteps = steps;
opt.silent = true;

out = spacarlight(nodes,elements,nprops,eprops,opt);


% addpath("C:\dvogureckiy99\2_Education\ITMO\DMS\SPACAR")
% % clear; clc;
% %addpath('spacar')
% E      = 2.95*1e9;    %Young's modulus 2.95*1у9
% nu = 0.39; % Puasson coefficient
% G = 0.5*E/(1+nu); %shear modulus  
% % rho    = 1270;        %density
% % yield_stress = 76.6; 

% % L      = 20e-3;       %length [20 60]1e3 with step 1e-4                 all 400    
% % H      = 10e-3;       %width  [10 50]1e3 with step 1e-4                 all 400
% % W      = 0.0005;      %thickness from 5*4e-4 with step [5 25]4e-4    all 20
% % Stroke = deg2rad(70); %full stroke                                          70 always
% % steps = 3;           %loadsteps 

% N = 11;
% nodes = [linspace(0,L,N)',zeros(N,1),zeros(N,1)];
% n = size(nodes,1);
% elements = [];
% for i=1:n-1
%     elements = [elements;i,i+1];
% end
% nprops(1).fix   = true; 
% F=30;
% l_Fext = 0.01;
% if exist('F','var')
% nprops(int32((l_Fext/L)*N)).force = [0,F,0];
% end
% % nprops(n).rot_z = Stroke;

% %flexures
% eprops(1).elems = linspace(1,n-1,n-1);
% eprops(1).flex = [3 4];
% % eprops(1).nbeams = 1;
% eprops(1).emod = E;
% eprops(1).smod = G;
% eprops(1).dens = double(rho);
% eprops(1).cshape ='rect';
% eprops(1).dim = [H W];
% eprops(1).orien = [0 0 1];
% eprops(1).color = 'blue';
% eprops(1).cw = true;

% opt.filename = 'simple_beam';
% opt.loadsteps = steps;
% opt.silent = true;

% out = spacarlight(nodes,elements,nprops,eprops,opt);

x_all_nodes = [];
y_all_nodes = [];
for i=1:n
    x_all_nodes = [x_all_nodes;out.node(i).p(1,steps)];
    y_all_nodes = [y_all_nodes;out.node(i).p(2,steps)];
end
phi_end = rad2deg(out.node(n).r_eulzyx(1,end));


