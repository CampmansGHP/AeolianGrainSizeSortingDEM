function [FLAGS,NUMSET,PARAMS,r,v,omega,d] = input_function
    Shield = 0.05;
    PARAMS.V_bridge_star=0; % 0.001;

    FLAGS.loadparticles = 2;    % 0: generate new results
                                % 1: load from INPUT dir
                                % 2: load from RUNS dir
                                % 3: generate from sediment sample
    load_name =  'example_generating_bed'; 
    %load_particle_sample_name = 'sediment_samples_P1500_d0_0.25_std_phi_0.5';
    
    % properties required for timestep settings.
    PARAMS.g = 9.81; % gravitational acceleration;
    % particle properties.
    PARAMS.rho_p = 2700;                                    % particle density;
    PARAMS.d0    = 0.00025;                                 % particle diameter;
    PARAMS.phi_std = 0.5;
    PARAMS.V0    = 4/3*pi*(PARAMS.d0/2)^3;                  % reference volume of particle;
    PARAMS.m0    = PARAMS.V0*PARAMS.rho_p;                  % reference mass of particle;
    PARAMS.k     = 5000*PARAMS.m0*PARAMS.g/PARAMS.d0;        % spring stiffness;  [Duran etal 2012, beneith eq (28)]
    PARAMS.nu_n  = 4.7*PARAMS.m0*sqrt(PARAMS.g/PARAMS.d0);    % damping coeficient normal direction; [Duran etal 2012, beneith eq (28)]
%     PARAMS.tc = pi/sqrt(2*PARAMS.k/PARAMS.m0 - (PARAMS.nu_n/PARAMS.m0)^2); % contact time; [Duran etal 2012, beneith eq (28), Luding (1998)]
    
    % time numerical settings; 
    NUMSET.t0    = 0;
    PARAMS.Tmax  = 2; 
    % (part of these settings are based on the smallest grain) (below the sediment grain size section.)
      
    PARAMS.P     = 40;%0;%150;        % number of particles; (including fixed particles at bed)
    FLAGS.particles_fixed_at_bed = 1; % true (1) or false (0)
    FLAGS.flow = 1;                   % true (1) or false (0)
    
    PARAMS.Ldx = 10;%0;%20;
    PARAMS.Lz  = 0.2;           % Height of domain;
    
    
    if FLAGS.loadparticles==1
        %% load initial conditions;
%         load(['../',load_name,'/r_vec.mat'],'r_vec');
%         load(['../',load_name,'/omega_vec.mat'],'omega_vec');
%         load(['../',load_name,'/v_vec.mat'],'v_vec');
        load('r_vec_low_tres.mat','r_vec_low_tres');
%         load('omega_vec_low_tres.mat','omega_vec_low_tres');
        load('v_vec_low_tres.mat','v_vec_low_tres');
        load('d.mat','d');
        r=r_vec_low_tres(:,:,end-1);
        v=v_vec_low_tres(:,:,end-1);
        [~,P]=size(r);
        omega=zeros(3,P);
%         omega(3,:)=omega_vec(:,end-1);
        clear r_vec v_vec omega_vec;
    elseif FLAGS.loadparticles==2
        %% load initial conditions;
%         load(['../../RUNS/',load_name,'/r_vec.mat'],'r_vec');
%         load(['../../RUNS/',load_name,'/omega_vec.mat'],'omega_vec');
%         load(['../../RUNS/',load_name,'/v_vec.mat'],'v_vec');
        load(['../../RUNS/',load_name,'/r.mat'],'r');
%         load(['../../RUNS/',load_name,'/omega.mat'],'omega');
        load(['../../RUNS/',load_name,'/v.mat'],'v');

        load(['../../RUNS/',load_name,'/d.mat'],'d');
%         r=r_vec(:,:,end);
%         v=v_vec(:,:,end);
        [~,P]=size(r);
        omega=zeros(3,P);
%         omega(3,:)=omega(:,end);
        clear r_vec v_vec omega_vec;
    elseif FLAGS.loadparticles==0
        %% create initial conditions
        
        phi0 = d2phi(PARAMS.d0);
        phi  = phi0+PARAMS.phi_std*randn(1,PARAMS.P);
        d = phi2d(phi);
        if FLAGS.particles_fixed_at_bed
            d(1:PARAMS.Ldx)=PARAMS.d0;
        end            
%         d=ones(1,PARAMS.P)*PARAMS.d0;
%         if FLAGS.particles_fixed_at_bed 
%             d(PARAMS.Ldx+1:end) = abs(d(PARAMS.Ldx+1:end)  + randn(1,PARAMS.P-PARAMS.Ldx)*0.1*PARAMS.d0);
%         else
%             d = abs(d  + randn(1,PARAMS.P)*0.1*PARAMS.d0);
%         end
        
%         r_fixed = zeros(3,PARAMS.Ldx);
%         r_fixed(1,:)=0.5*PARAMS.d0:PARAMS.d0:PARAMS.Ldx*PARAMS.d0;
        current_dir = cd; cd ../../
        %                               (P       ,Lx                  ,Ly         ,d)
        [r,v,omega] = initiate_particles(PARAMS.P,PARAMS.Ldx*PARAMS.d0,PARAMS.Lz/3,d);
        cd(current_dir);
        if FLAGS.particles_fixed_at_bed 
            r(1,1:PARAMS.Ldx) = 0.5*PARAMS.d0:PARAMS.d0:PARAMS.Ldx*PARAMS.d0;
            r(2,1:PARAMS.Ldx) = 0;
            v(:,1:PARAMS.Ldx) = 0;
            omega(:,1:PARAMS.Ldx) = 0;
        end
    elseif FLAGS.loadparticles==3
        %% generate initial positions of particles for a specified sample.
        load(['../../sediment_samples/',load_particle_sample_name,'.mat'],'sediment_size_data');
        if PARAMS.d0~=sediment_size_data.d0 || PARAMS.P~=numel(sediment_size_data.d)+PARAMS.Ldx
            disp('warning! The input d0 size or the number of particles does not match the sediment sample that is being loaded.');
            pause;
        end
        d = ones(1,PARAMS.P)*sediment_size_data.d0;
        d(PARAMS.Ldx+1:end)=sediment_size_data.d;
        current_dir = cd; cd ../../
        %                               (P       ,Lx                  ,Ly         ,d)
        [r,v,omega] = initiate_particles(PARAMS.P,PARAMS.Ldx*PARAMS.d0,PARAMS.Lz/3,d);
        cd(current_dir);
        if FLAGS.particles_fixed_at_bed 
            r(1,1:PARAMS.Ldx) = 0.5*PARAMS.d0:PARAMS.d0:PARAMS.Ldx*PARAMS.d0;
            r(2,1:PARAMS.Ldx) = 0;
            v(:,1:PARAMS.Ldx) = 0;
            omega(:,1:PARAMS.Ldx) = 0;
        end
%         whos r v omega d
%         pause;
    end
    
    m_smallest = PARAMS.rho_p*4/3*pi*(min(d)/2)^3;
    PARAMS.tc = pi/sqrt(2*PARAMS.k/m_smallest - (PARAMS.nu_n/m_smallest)^2); % contact time; [Duran etal 2012, beneith eq (28), Luding (1998)]
    NUMSET.dt    = 0.25*0.2*PARAMS.tc; 
    NUMSET.It    = ceil(PARAMS.Tmax/NUMSET.dt);  
    
    
%     d(:)=PARAMS.d0;
%     r(:,1) = [3*PARAMS.d0 1*PARAMS.d0 0];
%     r(:,2) = [4.1*PARAMS.d0 3*PARAMS.d0 0];
%     v(:,1) = [0.01 0 0];
%     v(:,2) = [-0.1 0 0];
    
    
%     PARAMS.ywall = 0.0;
    
    FLAGS.realtime_animation = 1;               % plot during computation.
    FLAGS.saveresults = 1; % [ 0 1 2 ] 0: dont save, 1: save processed data, 2: save each particle position.
                            % go for option 1 over option 2: saving the
                            % full partilce trajectories quickly generates
                            % enourmous save files.
    FLAGS.showomega3 = 0;
    FLAGS.plot_numbers = 0;
    FLAGS.plot_numbers_t0 = 0;
    FLAGS.showEnergy = 1;
    FLAGS.showomega3_end_of_sim = 0;
    FLAGS.plot_trajectories = 0;
    
    FLAGS.display_tauP_over_tau_star = 0;
    
    % store_name = ['tmp_P',num2str(P),'_result'];%'P40_d0p01_fluid';
%     store_name = ['tmp_P',num2str(P),'_result_ustar3'];%'P40_d0p01_fluid';

    

    % fluid parameters;
    PARAMS.rho_f  = 1.225;  % [kg/m^3]
    PARAMS.nu_f   = 1.4e-5; % [m^2/s];
    PARAMS.Cd_inf = 0.5;

    
    
    % particle parameters;
%     PARAMS.Emod  = 31.14e6;  % GPa (not used?)
    PARAMS.q      = 2/5;
    PARAMS.u_star = sqrt(Shield*(PARAMS.rho_p-PARAMS.rho_f)*PARAMS.g*PARAMS.d0/PARAMS.rho_f);
    PARAMS.I0     = PARAMS.q*PARAMS.m0*(PARAMS.d0/2).^2;     % reference rotational inertia;
%     PARAMS.Emod  = 31.14e6;  % GPa
    NUMSET.dz_min = 0.5*min(0.2*PARAMS.d0,PARAMS.nu_f/PARAMS.u_star);
    NUMSET.Iz     = 1000; % flow resolution;

    PARAMS.nu_t  = 0.01*PARAMS.m0*sqrt(PARAMS.g/PARAMS.d0); % damping coeficient tangential direction;
    PARAMS.mu    = 0.5;               % tangential force is not allowed to exceed mu times the normal force; ft = min(ft,mu*fn);
    PARAMS.mu_r  = 0;% 1e3*PARAMS.I0;%1e-2;%1e-6;%0.1;%0.1;%1;
    PARAMS.omega_max_r =inf;%2*pi;%2*pi;%inf;% 2*pi;
    
    
    if FLAGS.flow == 0
        PARAMS.Lz_plot    = PARAMS.P/PARAMS.Ldx*PARAMS.d0*2;
    else
        PARAMS.Lz_plot    = PARAMS.P/PARAMS.Ldx*PARAMS.d0*4;
    end
end

function phi = d2phi(d)
    phi = -log(d*1000)/log(2);
end

function d = phi2d(phi)
    d = 2.^(-phi)/1000;
end