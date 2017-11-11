%% Plane_Wave_FE_Code
% ======================================================================= %
% Dimitri Krattiger
% 10-12-2012
close all
clear; clc; tic; format compact
% profile clear
% profile on
% warning off all

% Server parameters
addpath(genpath('libraries'))

%% Check for & Create Save directories
% ======================================================================= %

save_results = false;
load_results = false;

% if save_data folder does not exist, make it
if save_results
    if ~exist('save_data','dir')
        mkdir('save_data');
    end
    if ~exist('save_data/models','dir')
        mkdir('save_data/models');
    end
    if ~exist('save_data/solutions','dir')
        mkdir('save_data/solutions');
    end
    if ~exist('save_data/BMS_solutions','dir')
        mkdir('save_data/BMS_solutions');
    end
    if ~exist('figures','dir')
        mkdir('figures');
    end
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

% Thumbnail Image Generator (Note that the videos and figures folder
% contains a much more capable thumbnail plotting code)
% mrl = 1; n = 8; cir_dist = 0;% (1960 DOF) 


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
        
        % Full FE model w(k) Dispersion Solution
        tic
        [omega_full,PHI_full,t_kloop_full] = dispersion_solver_w_k(kappa,K_free,M_free,dof_sets,R,n_curves);
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
options_BMS.InteriorMethod  = 'CB';
options_BMS.BoundaryMethod  = 'hybrid';
% options_BMS.w_i             = 2*w_cut_max;
% options_BMS.w_b             = 1*w_cut_max;
options_BMS.n_FI             = 200;
options_BMS.n_CC             = 22;
% options_BMS.n_CC             = 12;
options_BMS.verbose         = true;
options_BMS.plots           = true;

% perform BMS reduction 
[K_BMS,M_BMS,dof_sets_BMS,info_BMS,T_BMS] = BMS(K_free,M_free,coordinates,R,options_BMS);

% compute dispersion
[w_BMS,PHI_BMS,t_kloop_BMS] = dispersion_solver_w_k(kappa,K_BMS,M_BMS,dof_sets_BMS,R,n_curves);

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

% pause
%% Perform BMS Plus Reduction
% =================================================================== %



% interior reduction parameters
clear options_BMSpl
options_BMSpl.InteriorMethod       = 'CB+';
options_BMSpl.BoundaryMethod        = 'exact';
options_BMSpl.n_FI                  = 30;
options_BMSpl.n_CC                  = 12;

% AMLS+ interior reduction parameters
% % % %
% clear options_BMSpl
% options_BMSpl.InteriorMethod        = 'AMLS+';
% options_BMSpl.BoundaryMethod        = 'none';
% % options_BMSpl.n_FI                  = 30;
% % options_BMSpl.n_CC                  = 12;
% options_BMSpl.w_i                   = max(w_BMS(:));



% options_BMSpl.BoundaryMethod        = 'hybrid';
% options_BMSpl.n_CC                  = 25;


% options_BMSpl.BoundaryMethod        = 'exact';
% options_BMSpl.n_CC                  = 14;

% perform BMS reduction 
% [K_BMSpl,M_BMSpl,dof_sets_BMS,t_up_front_plus,T_BMS] = BMS_plus(K_free,M_free,coordinates,w_cut*0.5,R,n_FI,n_LI);
[K_BMSpl,M_BMSpl,dof_sets_BMSpl,info_BMSpl,T_BMSpl] = BMS(K_free,M_free,coordinates,R,options_BMSpl);

% solve for BMS dispersion
[w_BMSpl,PHI_BMSpl,t_kloop_BMSpl] = dispersion_solver_w_k(kappa,K_BMSpl,M_BMSpl,dof_sets_BMSpl,R,n_curves);

% convert to Hz
f_BMS_plus = w_BMSpl/(2*pi);

% evaluate and display frequency errors
tol = 1e-5;
e_freq2 = 100*abs(f_BMS_plus-f_full)./f_full;
e_freq_max2 = max(max(e_freq2(f_full>tol*max(max(f_full)))))
t_BMS_plus = sum(t_kloop_BMSpl) + info_BMSpl.t_up_front

t_full/sum(t_kloop_BMSpl)
t_full/t_BMS_plus


figure(33);clf
plot3(kappa_plot,f_BMS_plus*2*pi*norm_fac,e_freq2)

% plot dispersion
figure(4);clf('reset')
dispersion_plot(kappa_plot,{f_full*2*pi*norm_fac;f_BMS_plus*2*pi*norm_fac},...
    {'\Gamma','X','M','\Gamma'},{'Full Solution','BMS+ Solution'});
drawnow
% pause

%% plot Transformation modes
% ======================================================================= %
plot_transform = false;
if plot_transform


    % choose mode to plot
    i_phi_plot = 27;
    i_phi_plot = 34;
%     i_phi_plot = 1;
    
%     % fun modes to use for residual
%     i_phi_plot = 23;
%     i_phi_plot = 30;
%     i_phi_plot = 33;
%     i_phi_plot = 39;
%     i_phi_plot = 46;
    
    
    
    for i_phi_plot = 50:80
        
        % mode shape magnitude
        mag_PHI = 0.5;
        PHI_plot = T_BMS(:,i_phi_plot);
        PHI_plot = (Lx*mag_PHI)*PHI_plot/norm(PHI_plot);

        % Add mode shapes to original coordinates
        coordinates_phi = coordinates;
        coordinates_phi(:,1) = coordinates_phi(:,1) + real(PHI_plot(1:2:end));
        coordinates_phi(:,2) = coordinates_phi(:,2) + real(PHI_plot(2:2:end));

        figure(22);clf
        plot_FEM_model(coordinates_phi,patchfaces,C,fedges);
        axis equal
        drawnow
        i_phi_plot
        pause
    end
end

%% Evaluate BMS Plus Reduction
% ======================================================================= %
eval_BMS = true;
if eval_BMS
    
    % Interior Methods
    InteriorMethods = {'CB','CB+'};
    linestyles = {'-','--',':','-.'};
    colors = get(groot,'factoryAxesColorOrder');
    
    % full run values
    run_select = 'full';
    switch run_select
        case 'full' % full run values 
            
            % number of fixed-interface modes
            n_FI_list = ceil(logspace(log10(5),log10(500),10));
            n_FI_cell = num2cell(n_FI_list);
            
            % number of boundary CC modes
            CC_multiplier = 0.25;
            n_CC_list = ceil(n_FI_list*CC_multiplier);
            n_CC_cell = num2cell(n_CC_list);
            
            [runParam1(1:length(n_FI_list)).n_FI] = deal(n_FI_cell{:});
            [runParam1(1:length(n_CC_list)).n_CC] = deal(n_CC_cell{:});
            
            % interior reduction types  
            [runParam2(1:2).InteriorMethod] = deal('CB','CB+');            
            
            % third run parameter placeholder
            [runParam3(1:2).BoundaryMethod] = deal('none','exact');
            
            % set fixed parameters
            clear options_BMS;            
            solution_type = 'bandstructure';
            
        case 'small' % small test values
     
            % number of fixed-interface modes
            n_FI_list = ceil(logspace(log10(5),log10(100),5));
            n_FI_cell = num2cell(n_FI_list);
            
            % number of boundary CC modes
            CC_multiplier = 0.25;
            n_CC_list = ceil(logspace(log10(5),log10(100),5)*CC_multiplier);
            n_CC_cell = num2cell(n_CC_list);
            
            [runParam1(1:length(n_FI_list)).n_FI] = deal(n_FI_cell{:});
            [runParam1(1:length(n_CC_list)).n_CC] = deal(n_CC_cell{:});
            
            % interior reduction types  
            [runParam2(1:2).InteriorMethod] = deal('CB','CB+');            
            
            % third run parameter placeholder
            [runParam3(1:2).BoundaryMethod] = deal('none','exact');
            
            % set fixed parameters
            clear options_BMS;            
            solution_type = 'bandstructure';

        case 'AMLS_vs_HCB'
            
            % cutoff frequencies
            w_i_list = max(f_full(:))*2*pi*logspace(log10(3),1,5);
            w_i_cell = num2cell(w_i_list);
            [runParam1(1:5).w_i] = deal(w_i_cell{:});
            
            % interior reduction types  
            [runParam2(1:2).InteriorMethod] = deal('AMLS+','CB+');
            
            
            % third run parameter 
            % placeholder if not needed:   runParam3(1).dummy = [];
            [runParam3(1).BoundaryMethod] = deal('none');
            
            % set fixed parameters
            clear options_BMS;
            options_BMS.BoundaryMethod = 'none';
            
            solution_type = 'free';   
            
    end
    
    
    % number of evaluation parameters (along different dimensions)
    np1 = length(runParam1);
    np2 = length(runParam2);
    np3 = length(runParam3);

    % frequencies to check in Dispersion calculation
    i_check = abs(f_full)>1e-1;

    % preallocate arrays
    n_FI_save               = nan(np1,np2,np3);
    n_LI_save               = nan(np1,np2,np3);
    n_dof_BMS_save          = nan(np1,np2,np3);
    e_f_maxs                = nan(np1,np2,np3);
    
    t_BMSs                  = zeros(np1,np2,np3);
    t_BMS_ks                = zeros(np1,np2,np3);
    t_kloop_BMSs            = zeros(np1,np2,np3,n_kap);
    t_up_fronts             = zeros(np1,np2,np3); 
    
    legendstrings_BMS       = cell(np1*np2*np3,1);
    
    for k = 1:np3
        for j = 1:np2
            for i = 1:np1
                
                % each runParam structure contains a set of fields that
                % change in unison. 
                fnames1 = fieldnames(runParam1);
                fnames2 = fieldnames(runParam2);
                fnames3 = fieldnames(runParam3);
                
                % Setup BMS options
                for ii = 1:length(fnames1)
                    options_BMS.(fnames1{ii}) = runParam1(i).(fnames1{ii});
                end
                
                for jj = 1:length(fnames2)
                    options_BMS.(fnames2{jj}) = runParam2(j).(fnames2{jj});
                end
                
                for kk = 1:length(fnames3)
                    options_BMS.(fnames3{kk}) = runParam3(k).(fnames3{kk});
                end

                % Define a Save path
                % All this just define a string and path for saving?!
                % sadly yes
                BMSsolutionpathstring = 'save_data/BMS_solutions/';

                solutiondescription = sprintf('%ikpts_%iBands',n_kap,n_curves);
                    
                interior_type = [options_BMS.InteriorMethod,'_'];
                if isfield(options_BMS,'n_FI')
                    interior_size = sprintf('%iFImodes_',options_BMS.n_FI);
                else
                    interior_size = sprintf('%4.3eHzFreqCutoff_',options_BMS.w_i);
                end
                    
                if strcmpi(options_BMS.BoundaryMethod,'none')
                    boundary_type = 'NoBoundReduction';
                    boundary_size = [];
                    boundary_orthog = [];
                else
                    boundary_type = [upper(options_BMS.BoundaryMethod(1)),lower(upper(options_BMS.BoundaryMethod(2:end))),'BoundReduct_'];
                    boundary_size = sprintf('%iInitialCCModes_',options_BMS.n_CC);
                    if isfield(options_BMS,'orthoTypeLIRExact')
                        boundary_orthog = ['BoundOrthog_',upper(options_BMS.orthoTypeLIRExact)];
                    else
                        boundary_orthog = 'BoundOrthogQR';
                    end
                end
                    
                BMSdescription = [interior_type,interior_size,boundary_type,boundary_size,...
                                  boundary_orthog];                              
                              
                BMSsolution_savestring = [BMSsolutionpathstring,modeldescription,...
                                '_',solutiondescription,'_',BMSdescription];
                
                            
                if exist([BMSsolution_savestring,'.mat'],'file')  && load_results
                    load([BMSsolution_savestring,'.mat'])
                else
                            
                    %[K_BMSpl,M_BMSpl,dof_sets_BMSpl,t_up_front_plus_save(i,j),T_BMS_plus] = BMS_plus(K_free,M_free,coordinates,R,opts_BMS);
                    [K_BMS2,M_BMS2,dof_sets_BMS2,info_BMS2,T_BMS2] = BMS(K_free,M_free,coordinates,R,options_BMS);
                    %[K_BMS,M_BMS,dof_sets_BMS,t_up_fronts(k,i,j)] = BMS(K_free,M_free,coordinates,R,options_BMS);
                    
                    % compute BMS dispersion
                    [w_BMS2,PHI_BMS2,t_kloop_BMS2] = ...
                     dispersion_solver_w_k(kappa,K_BMS2,M_BMS2,dof_sets_BMS2,R,n_curves);

                    % convert to Hz
                    f_BMS2 = w_BMS2/(2*pi);
                          
                    % save BMS results
                    if save_results
                        save([BMSsolution_savestring,'.mat'],'f_BMS2','PHI_BMS2','t_kloop_BMS2','info_BMS2')
                    end
                end
                
                % store results
                f_BMS_save{i,j,k} = f_BMS2;
                PHI_BMS_save{i,j,k} = PHI_BMS2;
                t_kloop_BMSs(i,j,k,:) = t_kloop_BMS2;
                t_up_fronts(i,j,k) = info_BMS2.t_up_front;
                
                % timing data
                t_BMSs(i,j,k) = sum(t_kloop_BMSs(i,j,k,:)) + t_up_fronts(i,j,k);
                t_BMS_ks(i,j,k) = sum(t_kloop_BMSs(i,j,k,:))/n_kap;

                % model dimensions
                n_FI_save(i,j,k) = info_BMS2.n_FI;
                n_dof_BMS_save(i,j,k) = size(PHI_BMS2,1);
                n_LI_save(i,j,k) = n_dof_BMS_save(i,j,k)-n_FI_save(i,j,k);

                % evaluate error and store maximum error
                e_f_BMS2 = zeros(size(f_full));
                e_f_BMS2(i_check) = 100*(f_BMS2(i_check)-f_full(i_check))./f_full(i_check);
                e_f_maxs(i,j,k) = max(max(abs(e_f_BMS2)));


                % create legend strings
                legendstrings_BMS{(k-1)*np2+j} = [options_BMS.InteriorMethod,' BoundaryReduction = ',options_BMS.BoundaryMethod];

                % plot intermediate results
                figure(12);clf
                for kk = 1:np3
                    for jj = 1:np2
                        h1 = semilogy(squeeze(t_BMSs(:,jj,kk))/t_full,squeeze(e_f_maxs(:,jj,kk)));hold on
%                         h1 = semilogy(squeeze(n_dof_BMS_save(:,jj,kk)),squeeze(e_f_maxs(:,jj,kk)));hold on
                        
                        set(h1,'marker','.',...
                               'linestyle',linestyles{kk},...
                               'linewidth',2,...
                               'markersize',12,...
                               'color',colors(jj,:));
                        xlabel('Time Fraction, R');ylabel('Maximum Frequency Error (%)')
                    end
                end
                xlim([0,1])
                set(gca,'XMinorTick','off','YMinorTick','off')
                legendstrings_temp = legendstrings_BMS(1:((k-1)*np2+j));
                legend(legendstrings_temp{:},'location','northeastoutside')
                drawnow
                
                if save_results
                    savestring = sprintf('BMS Error Analysis %iDOF',n_dof_per);
                    saveas(gcf,['figures/',savestring,date_string])
                end

                % output some model-size info
                n_FI_save
                n_LI_save
            end
        end
    end
end

n_dof_BMS_save(:,1)./n_dof_BMS_save(:,2)
pause
%% plot mode
% ======================================================================= %

plot_mode = true;
if plot_mode

    % choose mode to plot
    i_phi_plot = 7;
    k_sel = 9;

    % mode shape magnitude
    mag_PHI = 0.07;
    PHI_plot = PHI_full(:,i_phi_plot,k_sel);
    [~,i_max] = max(imag(PHI_plot));
    PHI_plot = PHI_plot/PHI_plot(i_max);
    PHI_plot = PHI_plot/max(abs(PHI_plot));
    PHI_plot = (Lx*mag_PHI)*PHI_plot;
    
    % expand mode with periodicity transformation
    kvec = kappa(:,k_sel);
    lam  = exp(-1i*R*kvec);    
    lam(end+1:3) = 0;
        
    T_per = Periodic_Boundary_Conditions(dof_sets);
    T_per_k = T_per.s0 + T_per.s1*lam(1) + T_per.s2*lam(2) + T_per.s3*lam(3)...
            + T_per.s12*lam(1)*lam(2) + T_per.s23*lam(2)*lam(3) + T_per.s13*lam(1)*lam(3)...
            + T_per.s123*lam(1)*lam(2)*lam(3);
    PHI_plot = T_per_k*PHI_plot;
        
    % Add mode shapes to original coordinates
    coordinates_phi = coordinates;
    coordinates_phi(:,1) = coordinates_phi(:,1) + real(PHI_plot(1:2:end));
    coordinates_phi(:,2) = coordinates_phi(:,2) + real(PHI_plot(2:2:end));

    PHI_plot2 = BMS_Plus_Mode_Expansion(PHI_BMSpl(:,i_phi_plot,k_sel),...
        dof_sets_BMSpl,kvec,R,T_BMSpl,K_BMSpl,M_BMSpl);

    PHI_plot2 = PHI_plot2*(PHI_plot2'*PHI_plot);

    PHI_plot2 = (Lx*mag_PHI)*PHI_plot2/[diag(max(abs(PHI_plot2)))]';
    % PHI_plot2 = sum(PHI_plot2,2);
    coordinates_phi2 = coordinates;
    coordinates_phi2(:,1) = coordinates_phi2(:,1) + real(PHI_plot2(1:2:end));
    coordinates_phi2(:,2) = coordinates_phi2(:,2) + real(PHI_plot2(2:2:end));

    % figure(123);clf
    % plot(coordinates_phi(:,1),'r.-');hold on
    % plot(coordinates_phi2(:,1),'bo--')
    % drawnow

    figure(22);clf
    % subplot(1,2,1)
    plot_FEM_model(coordinates_phi,patchfaces,C,fedges);
    axis equal

    % subplot(1,2,2)
    [hpatch,hline] = plot_FEM_model(coordinates_phi2,patchfaces,C,fedges);
    axis equal
    delete(hpatch)
    set(hline,'linestyle','--','color','green')

    drawnow

    % e_phi = sin(acos((PHI_plot'*PHI_plot2)/(sqrt(PHI_plot'*PHI_plot)*sqrt(PHI_plot2'*PHI_plot2))))
    e_phi = 1-abs(PHI_plot'*PHI_plot2)/(sqrt(PHI_plot'*PHI_plot)*sqrt(PHI_plot2'*PHI_plot2))
    % e_phi = acos(abs(PHI_plot'*PHI_plot2)/(sqrt(PHI_plot'*PHI_plot)*sqrt(PHI_plot2'*PHI_plot2)))*pi/180
    % e_phi = norm(PHI_plot/norm(PHI_plot)-PHI_plot2/norm(PHI_plot2))
end

%% evaluate maximum mode error
% ======================================================================= %
eval_modes = true;
if eval_modes
    T_per3 = Periodic_Boundary_Conditions(dof_sets_BMSpl);
    T_per2 = Periodic_Boundary_Conditions(dof_sets_BMS);
    e_phi = zeros(n_curves,n_kap);
    e_phi_plus = zeros(n_curves,n_kap);
    for i = 1:n_curves
        for j = 1:n_kap
            [i,j]            
            
            % phase modification at boundaries due to wavevector
            kvec = kappa(:,j);
            lam  = exp(-1i*R*kvec);

            lam(end+1:3) = 0;
            
            T_per_k = T_per.s0 + T_per.s1*lam(1) + T_per.s2*lam(2) + T_per.s3*lam(3)...
                    + T_per.s12*lam(1)*lam(2) + T_per.s23*lam(2)*lam(3) + T_per.s13*lam(1)*lam(3)...
                    + T_per.s123*lam(1)*lam(2)*lam(3);
            T_per_k2 = T_per2.s0 + T_per2.s1*lam(1) + T_per2.s2*lam(2) + T_per2.s3*lam(3)...
                    + T_per2.s12*lam(1)*lam(2) + T_per2.s23*lam(2)*lam(3) + T_per2.s13*lam(1)*lam(3)...
                    + T_per2.s123*lam(1)*lam(2)*lam(3);
                
                
            % extract modes
            phi_1 = T_per_k*PHI_full(:,i,j);
            phi_2 = T_BMS*T_per_k2*PHI_BMS(:,i,j);
            phi_3 = BMS_Plus_Mode_Expansion(PHI_BMSpl(:,i,j),dof_sets_BMSpl,kappa(:,j),R,T_BMSpl,K_BMSpl,M_BMSpl);
            
            % check for degeneracy (don't count error for degenerate modes)
            f_close = (abs(f_full(i,j)-f_full(:,j))/abs(f_full(i,j)));
            f_close(i) = nan;
            if any(f_close<1e-2) | abs(f_full(i,j))<1e-2*max(abs(f_full(:)));
                e_phi(i,j) = 0;
                e_phi_plus(i,j) = 0;
            else
    %             e_phis(i,j) = 1-abs(phi_1'*phi_2)/(sqrt(phi_1'*phi_1)*sqrt(phi_2'*phi_2));
                e_phi(i,j) = 1-abs(phi_1'*phi_2)/(sqrt(phi_1'*phi_1)*sqrt(phi_2'*phi_2));
                e_phi_plus(i,j) = 1-abs(phi_1'*phi_3)/(sqrt(phi_1'*phi_1)*sqrt(phi_3'*phi_3));

            end
        end
    end

    e_phi_max = max(max(e_phi))
    e_phi_plus_max = max(max(e_phi_plus))
    figure(212);clf
    plot(f_full,log10(e_phi),'b.');hold on
    plot(f_full,log10(e_phi_plus),'r.')
    
    xlabel('frequency'),ylabel('mode error')
    
end