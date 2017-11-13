%% Plane_Wave_FE_Code
% ======================================================================= %
% Dimitri Krattiger
% 10-12-2012
close all
clear; clc; tic; format compact
% profile clear
% profile on
warning off all

% Server parameters
addpath(genpath('Global_Files'))
addpath('Global_Files')

%% Check for & Create Save directories
% ======================================================================= %

save_results = true;
load_results = true;

% if save_data folder does not exist, make it
% if save_results
%     if ~exist('save_data','dir')
%         mkdir('save_data');
%     end
%     if ~exist('save_data/models','dir')
%         mkdir('save_data/models');
%     end
%     if ~exist('save_data/solutions','dir')
%         mkdir('save_data/solutions');
%     end
%     if ~exist('save_data/BMS_solutions','dir')
%         mkdir('save_data/BMS_solutions');
%     end
% end

if save_results
    mkdir('save_data');
    mkdir('save_data/models');
    mkdir('save_data/solutions');
    mkdir('save_data/BMS_solutions');
    mkdir('figures');    
end

%% Solution Options
% ======================================================================= %
j = 1;
% number of dispersion branches to compute
n_curves = 10;

% % number of fixed interface modes
% n_FI = 100;

% % number of characteristic constraint modes
% n_b_phi = 10;

% mode comparison
mode_plot = [1:10];
n_mp = length(mode_plot);
% k_sel = 9;

% models to compute
full_soln = true;
BMS_soln = true;
BMS_LIR_soln = true;

% save results?
% save_results = true;
% load_results = false;
date_string = datestr(now,'dd-mmm-yyyy HH-MM');

%% Circular Inclusion model dimensions
% ======================================================================= %

% unit cell dimensions
n_cellx = 1;
n_celly = 1;
Lx = 0.05*n_cellx;
Ly = 0.05*n_celly;

% lattice vectors
r1 = [Lx;0];
r2 = [0;Ly];
r3 = [];

R = [r1,r2];

% inclusion geometry
r_ss = Lx*0.1;  % small square radius
r_bs = Lx*0.5;  % big square radius
r_c = 1*Lx*0.25;    % circle radius
r_l = 1*Lx*0.25;    % layer radius

% % distortion parameters cube
cir_dist = 0.5; % distortion magnitude of circle inclusion
n_petals = 8;
theta_offset = 0*pi/180;

% Circle to plot FI modes
% cir_dist = 0; %
% sqr_dist = 0; %0.05
% n_petals = 8;
% theta_offset = 0*pi/180;
% sqr_bow = 0;
% 
% mrl = 3; n=3;

%% Define Pattern of Materials in Model
% ======================================================================= %

% mesh size
% mrl = 1; n=2; % (40 DOF)
% mrl = 2; n=2; % (160 DOF)
% mrl = 1; n=8; % (? DOF)
% mrl = 3; n=2; % (360 DOF)
% mrl = 4; n=2; % (640 DOF)
% mrl = 5; n=2; % (1000 DOF)
% mrl = 6; n=2; % (1440 DOF)
% mrl = 7; n=2; % (1960 DOF)
% mrl = 8; n=2; % (2560 DOF)
% mrl = 9; n=2; % (3240 DOF)
% mrl = 10; n=2; % (4000 DOF)
% mrl = 11; n=2; % (4840 DOF)
% mrl = 12; n=2; % (5760 DOF)
% mrl = 13; n=2; % (6760 DOF)
% mrl = 14; n=2; % (7840 DOF)
% mrl = 15; n=2; % (9000 DOF)
% mrl = 16; n=2; % (5112 DOF)
% mrl = 17; n=2; % (11560 DOF)
% mrl = 18; n=2; % (12960 DOF)
% mrl = 19; n=2; % (14440 DOF)
% mrl = 20; n=2; % (16000 DOF)
% mrl = 21; n=2; % (17640 DOF)
% mrl = 22; n=2; % (19360 DOF)
% mrl = 23; n=2; % (21160 DOF)
% mrl = 24; n=2; % (23040 DOF)
% mrl = 25; n=2; % (25000 DOF)

% mrl = 4; n=3; % (2560 DOF)
% mrl = 8; n=3; % (10240 DOF)
mrl = 12; n=3; % (23040 DOF)
% mrl = 16; n=3; % (40960 DOF)

% mrl = 4; n=4; % (5760 DOF)
% mrl = 8; n=4;% (23040 DOF)

% define number of elements in mesh based on mesh-refinement level (mrl)
n_ele_ss = 2*mrl; % number of elements along edge of small internal square
n_ele_c = 1*mrl;  % number of radial elements between small square and circle
n_ele_layer = 0*mrl; % number of radial elements in layer padding the inclusion
n_ele_bs = 1*mrl;    % number of radial elements between inclusion and unit cell edge

%% Element Properties
% ======================================================================= %

% Fiber
rho1    = 2.713e3;
E1      = 70e9;
v1      = 0.34;
lam1    = E1*v1/((1+v1)*(1-2*v1));
mu1     = E1/(2*(1+v1));
D1      = [lam1+2*mu1, lam1, 0; lam1,lam1+2*    mu1, 0; 0,0, mu1];

% Matrix
rho2    = rho1/8;
E2      = E1/16;
v2      = 0.34;
lam2    = E2*v2/((1+v2)*(1-2*v2));
mu2     = E2/(2*(1+v2));
D2      = [lam2+2*mu2, lam2, 0; lam2,lam2+2*mu2, 0; 0,0, mu2];

% Void
rho3    = rho1/8;
E3      = E1/16;
v3      = 0.34;
lam3    = E2*v2/((1+v2)*(1-2*v2));
mu3     = E2/(2*(1+v2));
D3      = [lam2+2*mu2, lam2, 0; lam2,lam2+2*mu2, 0; 0,0, mu2];

% Averaged Material
D3 = (D1+D2)/2;
rho3 = (rho1+rho2)/2;

% Perturbed Averaged Material 1 
perturbation  = 0.0001;
D4 = D3*(1+perturbation);
rho4 = rho3*(1+perturbation);

% Perturbed Averaged Material 2 
D5 = D3*(1-perturbation);
rho5 = rho3*(1-perturbation);

% element properties;
Ds = {D1,D2,D3,D4,D5};
rhos = {rho1,rho2,rho3,rho4,rho5};

%% Create Wave vector
% ======================================================================= %

sym_pts = {'Gamma','X','M','Gamma'};
kap_ref = 5;
[kappa,kappa_plot] = wave_vector(sym_pts,kap_ref,R);
n_kap = size(kappa,2);

%% Create Mesh
% ======================================================================= %

[xlocs,ylocs,elenodes,pattern,patchfaces,C,fedges] = ...
     mesh_circular_inclusion(n,n_ele_ss,n_ele_c,n_ele_layer,...
         n_ele_bs,r_ss,r_bs,r_c,r_l,cir_dist,n_petals,theta_offset);

%% Find node sets for boundaries
% ======================================================================= %

% compute node sets for overall unit cell
[node_sets] = find_node_sets([xlocs,ylocs],[r1,r2]);
dof_sets = node2dof(node_sets,2);

% %% define element connectivity matrix
% % ======================================================================= %
% corners = [1,n,n^2,n^2-n+1];
% sides = [(2:n-1),(n+n:n:n^2-n),(n^2-1:-1:n^2-n+2),(n^2-2*n+1:-n:n+1)];
% interior = 1:n^2; interior([corners,sides]) = [];
% elenode_order = [corners,sides,interior];

%% Display unit cell geometry
% ======================================================================= %

coordinates = [xlocs,ylocs];

figure(2);clf
h_patch = plot_FEM_model(coordinates,patchfaces,C,fedges);
axis equal
title('Unit Cell')
axis off
% set(h_patch,'linewidth',4)
drawnow
% break
% pause

%% Form Transformation to enforce full model periodicity
% ======================================================================= %
T_per = Periodic_Boundary_Conditions(dof_sets);

n_nodes = length(xlocs);
n_dof = size(T_per.s0,1);
n_dof_per = size(T_per.s0,2);
n_elements = size(elenodes,1);

%% Form Master Mass and Stiffness
% ======================================================================= %

% define model save path
pathstring = 'save_data/models/';

% define element order string
orderstrings = {'1st','2nd','3rd','4th','5th',...
                '6th','7th','8th','9th','10th'};   
orderstring = orderstrings{n-1};

% full model description string
modeldescription = [sprintf('%iLobeInclusion_%iDOF_',n_petals,n_dof_per),...
                    orderstring,'OrderEles'];

% model_savestring = sprintf('data/%i_lobe_inclusion_model_%iDOF_%iNodesPerEle',...
%                 n_petals,n_dof_per,n^2);

model_savestring = [pathstring,modeldescription,'_FreeModel'];% form master mass and stiffness matrices

            
if exist([model_savestring,'.mat'],'file') && load_results
    load([model_savestring,'.mat'])
else

    % store mesh info in a structure
    Mesh_info.Ds            = Ds;
    Mesh_info.rhos          = rhos;
    Mesh_info.pattern_vec   = pattern;
    Mesh_info.n             = n;
    Mesh_info.n_dof         = n_dof;
    Mesh_info.elenodes      = elenodes;
    Mesh_info.xlocs         = xlocs;
    Mesh_info.ylocs         = ylocs;

%     [K_free,M_free] = ...
%         master_mass_stiffness_pw(kappa(:,1),Mesh_info);
    [K_free,M_free] = master_mass_stiffness_planeStrainQuad(Mesh_info);
    
    if save_results
        save(model_savestring,'K_free','M_free','coordinates')
    end
end

%% normalization factor
% ======================================================================= %
norm_fac = sqrt(rho1/E1)*Lx;


%% w(k) Dispersion Solution
% ======================================================================= %

if full_soln
    % if full dispersion results exist for specified model size, then load
    % them. Otherwise compute them.
%     full_disp_savestring = sprintf('data/%i_lobe_inclusion_full_dispersion_%iDOF_%iNodesPerEle',...
%                 n_petals,n_dof_per,n^2);


    solutionpathstring = 'save_data/solutions/';

    solutiondescription = [sprintf('%i',n_kap),sprintf('kpts_%iBands',n_curves)];

    solution_savestring = [solutionpathstring,modeldescription,'_',solutiondescription];   

    if exist([solution_savestring,'.mat'],'file')  && load_results
        load([solution_savestring,'.mat'])
    else
        
        
        % solution options
        clear options;
        options.n_curves = n_curves;
        
        % Full FE model w(k) Dispersion Solution
        tic
        [omega_full,PHI_full,t_kloop_full] = dispersion_solver_w_k(kappa,K_free,M_free,dof_sets,R,options);
        t_full = toc;
        f_full = omega_full/(2*pi);
        
        if save_results
            save(solution_savestring,'f_full','PHI_full','t_kloop_full')
        end
    end
    
    % full computation time;
    t_full = sum(t_kloop_full);
end


%% Perform BMS Reduction
% =================================================================== %

% w_cut_guess = sqrt(E2/rho2)*(pi/sqrt(sum(sum(R,2).^2)))*(n_curves/2);
% w_cut_max = max(max(f_full*2*pi));
% w_cut_max = 5/(norm_fac);

clear options_BMS                
options_BMS.InteriorMethod   = 'CB';
options_BMS.BoundaryMethod   = 'exact';
% options_BMS.w_i             = 2*w_cut_max;
% options_BMS.w_b             = 1*w_cut_max;
options_BMS.n_FI             = 30;
options_BMS.n_CC             = 12;
options_BMS.verbose          = true;
options_BMS.plots            = true;

% perform BMS reduction 
[K_BMS,M_BMS,dof_sets_BMS,info_BMS,T_BMS] = BMS(K_free,M_free,coordinates,R,options_BMS);

% solution options
clear options;
options.n_curves = n_curves;

% compute dispersion
[w_BMS,PHI_BMS,t_kloop_BMS] = dispersion_solver_w_k(kappa,K_BMS,M_BMS,dof_sets_BMS,R,options);

% convert to Hz
f_BMS = w_BMS/(2*pi);

% plot dispersion
figure(3);clf('reset')
dispersion_plot(kappa_plot,{f_full*2*pi*norm_fac;w_BMS*norm_fac},...
    {'\Gamma','X','M','\Gamma'},{'Full Solution','BMS Solution'})
drawnow

% evaluate and display frequency errors
tol = 1e-3;
e_freq = 100*(f_BMS-f_full)./f_full;
e_freq_max = max(max(e_freq(f_full>tol*max(max(f_full)))))
t_BMS = sum(t_kloop_BMS) + info_BMS.t_up_front

%% Perform BMS Plus Reduction
% =================================================================== %


clear options_BMSpl

% interior reduction parameters
options_BMS_pl.InteriorMethod       = 'CB+';
options_BMSpl.n_FI                  = 30;
options_BMSpl.n_CC                  = 12;
options_BMSpl.verbose               = true;
options_BMSpl.plots                 = true;
options_BMSpl.orthoTypeLIRExact     = 'svd';
options_BMS_pl.verbose              = true;

% perform BMS reduction 
% [K_BMSpl,M_BMSpl,dof_sets_BMS,t_up_front_plus,T_BMS] = BMS_plus(K_free,M_free,coordinates,w_cut*0.5,R,n_FI,n_LI);
[K_BMSpl,M_BMSpl,dof_sets_BMSpl,info_BMS_plus,T_BMS_plus] = BMS(K_free,M_free,coordinates,R,options_BMSpl);

% solve for BMS dispersion
clear options;
options.n_curves = n_curves;
[w_BMS_plus,PHI_BMS_plus,t_kloop_BMS_plus] = dispersion_solver_w_k(kappa,K_BMSpl,M_BMSpl,dof_sets_BMSpl,R,options);

% convert to Hz
f_BMS_plus = w_BMS_plus/(2*pi);

% evaluate and display frequency errors
tol = 1e-5;
e_freq2 = 100*abs(f_BMS_plus-f_full)./f_full;
e_freq_max2 = max(max(e_freq2(f_full>tol*max(max(f_full)))))
t_BMS_plus = sum(t_kloop_BMS_plus) + info_BMS_plus.t_up_front

% plot dispersion
figure(4);clf('reset')
dispersion_plot(kappa_plot,{f_full*2*pi*norm_fac;f_BMS_plus*2*pi*norm_fac},...
    {'\Gamma','X','M','\Gamma'},{'Full Solution','BMS+ Solution'});
drawnow

%% Plot w(k) dispersion
% =================================================================== %

norm_fac = sqrt(rho1/E1)*Lx;
figure(3);clf('reset')
dispersion_plot(kappa_plot,{f_full*2*pi*norm_fac;w_BMS*norm_fac;w_BMS_plus*norm_fac},...
    {'\Gamma','X','M','\Gamma'},{'Full Solution','BMS(CB)(LI)','BMS(CB+)(LI)'})

% w_cut_guess = sqrt(E2/rho2)*(pi/sqrt(sum(sum(R,2).^2)))*(n_curves/2);
if exist('f_full')
    w_cut_max = max(max(f_full*2*pi))*1;
else
    w_cut_max = sqrt(E2/rho2)*(pi/sqrt(sum(sum(R,2).^2)))*(n_curves/2);
end

%% Full k(w) Dispersion Solution
% ======================================================================= %

n_om = 200;
omega = linspace(0,w_cut_max,n_om);

% profile clear
% profile on
if full_soln
    
    solutionpathstring = 'save_data/solutions/';
    solutiondescription = [sprintf('%i',n_om),'wpts'];
    solution_savestring = [solutionpathstring,modeldescription,'_',solutiondescription]; 

    if exist([solution_savestring,'.mat'],'file')  && load_results
        load([solution_savestring,'.mat'])
    else
        
        % Full FE model k(w) Dispersion Solution
        tic
        
        clear options_k_w
        options_k_w.n_curves = n_curves;
        
        [kappa_full,phi,t_wloop_full] = ...
            dispersion_solver_k_w(omega,K_free,[],M_free,dof_sets,R,options_k_w);
        t_full_k = toc;
        if save_results
            save(solution_savestring,'kappa_full','t_wloop_full')
        end
    end
end

%% Compute BMS dispersion using k-w method
% ======================================================================= %
% omega = linspace(0,w_cut_max,500);

% profile clear
% n_curves = [];


clear options_k_w
options_k_w.n_curves = n_curves;

[kappa_BMS,~,t_wloop_BMS] = dispersion_solver_k_w(omega,K_BMS,[],M_BMS,...
    dof_sets_BMS,R,options_k_w);

sum(t_wloop_BMS)


%% Compute BMS+ dispersion using k-w method
% ======================================================================= %

n_om = 200;
clear options_k_w
options_k_w.n_curves = n_curves;
omega_plus = linspace(0,w_cut_max,n_om);

[kappa_BMS_plus,PHI_BMS_plus,t_wloop_BMS_plus] = dispersion_solver_k_w(omega_plus,K_BMSpl,[],M_BMSpl,dof_sets_BMSpl,R,options_k_w);

% timing results
t_BMS_plus_k_w = sum(t_wloop_BMS_plus) + info_BMS_plus.t_up_front
t_full_k_w = sum(t_wloop_full)
t_full_k_w/t_BMS_plus_k_w

%% plot BMS k(w) solution
% ======================================================================= %
curve_plot = 1:8;

options_plot.Markers = {'o','v'};
options_plot.ThreeD = true;
options_plot.LineStyles = {'none','none'};
options_plot.Colors = {'k','g'};
om_plot = 1:length(omega);
omega_plot = omega(om_plot)*norm_fac;
kappas_plot = {kappa_full(curve_plot,om_plot),...
               kappa_BMS_plus(curve_plot,om_plot)};
figure(4);clf
[h,legendvec] = dispersion_plot_k_w(omega_plot,kappas_plot,options_plot);


% add in w(k) dispersion
hold on
k_plot = 1:((n_kap-1)/3+1);
% k_plot = 1:n_kap;
subplot(1,3,2)
h2 = plot(kappa_plot(k_plot),f_full(:,k_plot)*2*pi*norm_fac,'k-')

legend([legendvec,h2(1)],'BMS_{HCB+}','full','full_k(w)')

figure(5);clf;hold on;view(3)
h1 = plot3(real(kappa_BMS(curve_plot,:)),imag(kappa_BMS(curve_plot,:)),omega/(2*pi),'g.');hold on
h2 = plot3(real(kappa_BMS_plus(curve_plot,:)),imag(kappa_BMS_plus(curve_plot,:)),omega_plus/(2*pi),'bx');hold on
h3 = plot3(real(kappa_full(curve_plot,:)),imag(kappa_full(curve_plot,:)),omega/(2*pi),'ko');hold on
legend([h1(1),h2(1),h3(1)],'BMS','BMS+','Full')

%%
kappa_temp = kappa_BMS_plus(curve_plot,:);
max_real = max(abs(real(kappa_temp(:))));
max_imag = max(abs(imag(kappa_temp(:))));
i_real = abs(imag(kappa_temp))<max_imag*1e-3;
i_imag = abs(real(kappa_temp))<max_real*1e-3;
i_comp = ~i_real & ~i_imag;

kappa_real = kappa_temp;kappa_real(~i_real) = nan;
kappa_imag = kappa_temp;kappa_imag(~i_imag) = nan;
kappa_comp = kappa_temp;kappa_comp(~i_comp) = nan;

figure(99);clf
% plot3(abs(real(kappa_temp)),abs(imag(kappa_temp)),omega/(2*pi),'k.');hold on
plot3(abs(real(kappa_real)),-abs(imag(kappa_real)),omega_plus/(2*pi),'r.-');hold on
plot3(abs(real(kappa_imag)),-abs(imag(kappa_imag)),omega_plus/(2*pi),'g.-');hold on
plot3(abs(real(kappa_comp)),-abs(imag(kappa_comp)),omega_plus/(2*pi),'b.-');hold on
xlim([0,max_real*1.001])
ylim([-2*max_real,0])
set(gca,'xtick',[0,max_real*1.001])
set(gca,'ytick',max_real*[-2,-1,0])
grid on
xlabel('real')
ylabel('imaginary')
zlabel('frequency (Hz)')
set(gca, 'DataAspectRatio', [repmat(min(diff(get(gca, 'XLim')), diff(get(gca, 'YLim'))), [1 2]) diff(get(gca, 'ZLim'))])

figure(199);clf
plot(abs(real(kappa_temp)),omega_plus/(2*pi),'k.-');hold on
plot(-abs(imag(kappa_temp)),omega_plus/(2*pi),'k.-');hold on
xlim([-2*max_real,max_real])
% plot3(abs(real(kappa_comp)),-abs(imag(kappa_comp)),omega_plus/(2*pi),'b.-');hold on


%% evaluate BMS error
% ========================================================================
% %

kappa_BMS_squeeze = sort(abs(real(kappa_BMS(curve_plot,:))) + 1i*abs(imag(kappa_BMS(curve_plot,:))),1);
kappa_BMS_plus_squeeze = sort(abs(real(kappa_BMS_plus(curve_plot,:))) + 1i*abs(imag(kappa_BMS_plus(curve_plot,:))),1);
kappa_full_squeeze = sort(abs(real(kappa_full(curve_plot,:))) + 1i*abs(imag(kappa_full(curve_plot,:))),1);

for i = 1:n_om
    Cmat = abs(kron(kappa_BMS_plus_squeeze(:,i),ones(1,length(curve_plot)))-...
               kron(kappa_full_squeeze(:,i),ones(1,length(curve_plot))).') ; % cost matrix
        [i_sort] = match_points(Cmat); % dynamic programming matching
        kappa_BMS_plus_squeeze(:,i) = kappa_BMS_plus_squeeze(i_sort,i);
end
e_kap_BMS_plus = 100*max(max(abs(kappa_BMS_plus_squeeze - kappa_full_squeeze)))/(pi/Lx)
