function particle_model(FLAGS,PARAMS,NUMSET,r,v,omega,d,store_name)
    load('default_colors.mat','color');

    Nsubplots = FLAGS.showomega3+FLAGS.flow+1;
    Nflow = max(2,Nsubplots);
    Nomega = min(2,Nsubplots);
    
    d0      = PARAMS.d0;
    P       = PARAMS.P;
    g       = PARAMS.g;
    rho_f   = PARAMS.rho_f;
    rho_p   = PARAMS.rho_p;
    u_star  = PARAMS.u_star;
    q       = PARAMS.q;
    Ldx     = PARAMS.Ldx;
    nu_f    = PARAMS.nu_f;
    Cd_inf  = PARAMS.Cd_inf;
    
    Iz      = NUMSET.Iz;
    It      = NUMSET.It;
    t       = NUMSET.t0;
    dt      = NUMSET.dt;
    Tmax    = PARAMS.Tmax;
    
    V = 4/3*pi*(d/2).^3;
    m = V*rho_p;
    I = PARAMS.q*m.*(d/2).^2;
    
    if FLAGS.flow
        ag = [0;-PARAMS.g*(1-rho_f/rho_p);0];
    else
        ag = [0;-PARAMS.g;0];
    end

    Shields = rho_f*u_star^2/((rho_p-rho_f)*PARAMS.g*d0);
    disp(['Shields = ',num2str(Shields)]);

    Lx = d0*Ldx;
    % Lx = 0.025;     % Length of domain expressed in number of particles;
    Lz      = PARAMS.Lz;
    Lz_plot = PARAMS.Lz_plot;

%     A = Lx*d0;  
    
    k    = PARAMS.k;                    % spring stiffness;
    nu_n = PARAMS.nu_n;                 % damping coeficient normal direction;
    nu_t = PARAMS.nu_t;                 % damping coeficient tangential direction;
    mu   = PARAMS.mu;                   % tangential force is not allowed to exceed mu times the normal force; ft = min(ft,mu*fn);
    mu_r = PARAMS.mu_r;                 %
    omega_max_r = PARAMS.omega_max_r;   % 
    

   

%     r=r_vec(:,:,end);
%     v=v_vec(:,:,end);
%     [~,P]=size(r);
%     omega=zeros(3,P);
%     omega(3,:)=omega_vec(:,end);
%     clear r_vec v_vec omega_vec;
    % 
    % r(2,15)=Lz*0.5;
    % v(1,15)=0.1;

    
    %% plot t=0;
    figpart=figure('units','centimeters','position',[13 2 30 16]);
    subplot(1,Nsubplots,1);
    for p=1:P
                    myplot(r(:,p),d(p),'color',color(:,mod(p-1,7)+1));  % only these get a legend entry;
                    if p==1
                        hold on;
                    end
                    if FLAGS.plot_numbers || FLAGS.plot_numbers_t0; text(r(1,p),r(2,p),num2str(p)); end
                    if r(1,p)-Lx+d(p)/2>0
                        myplot(r(:,p)-[Lx;0;0],d(p),'color',color(:,mod(p-1,7)+1));
                    end
                    if r(1,p)+Lx-d(p)/2<Lx
                        myplot(r(:,p)+[Lx;0;0],d(p),'color',color(:,mod(p-1,7)+1));
                    end       
    end
    view(0,90);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('befor');
    axis equal;
    xlim([0 Lx]);
    ylim([0 Lz_plot]);
    drawnow;

    %% the sorting algorith seems to be very efficient in itself.
    % tic1=tic;
    if isfield(PARAMS,'V_bridge_star')
        V_bridge_star=PARAMS.V_bridge_star;
    else
        V_bridge_star = 0;
    end
    if V_bridge_star>0
        theta_contact_angle = 0;
        a_crit_star = (1+theta_contact_angle/2)*(V_bridge_star^(1/3)+V_bridge_star^(2/3)/10);
        a_crit = a_crit_star*d0;
    else
        a_crit = 0;
    end
    list = initiate_list(r,d+2*a_crit,Lx);
    % t1 = toc(tic1);
    % r(1,3) =0.8;

    % figure(2);
    % for p=1:P
    %                 myplot(r(:,p)+[Lx;0;0],d(p),'color',color(:,mod(p-1,7)+1));
    %                 if p==1
    %                     hold on;
    %                 end
    %     graph(p) =  myplot(r(:,p),d(p),'color',color(:,mod(p-1,7)+1));  % only these get a legend entry;
    %                 text(r(1,p),r(2,p),num2str(p));
    %                 myplot(r(:,p)-[Lx;0;0],d(p),'color',color(:,mod(p-1,7)+1));
    % end
    % view(0,90);
    % xlabel('x');
    % ylabel('y');
    % zlabel('z');
    % title('after');
    % axis equal;
    % xlim([0 Lx]);
    % ylim([0 Lz]);
    % drawnow;

    % tic2=tic;
    list = update_list(list,r,d+2*a_crit,Lx);
    % t2 = toc(tic2)

    %% Bounding Box Collisions

    % list
    % plot_bounding_boxes(list,color,Lx,Lz);

    % tic3=tic;
    % [delta,N] = collision_detection(list,d,r,Lx);
    % t3 = toc(tic3)


    tau_p_instant=zeros(Iz-1,1);
    if FLAGS.flow
        phi=zeros(Iz-1,1);
        [u,z,l,zph]=length_scale_model_function_particle_interaction(rho_f,nu_f,Lz,Iz,u_star,tau_p_instant,PARAMS.d0,phi); % initial profile;
        u_steady = u;
        U_plot_scale = max(u);
        if isfield(FLAGS,'flow_initiation')
            u = u.*(1-exp(-z/FLAGS.flow_initiation))';
    %         u=u*(FLAGS.flow_initiation+eps);
            Istartup_flow = 0;
        else
            Istartup_flow = 300;
        end
        f_fluid = calc_fluid_forces(u,z,r,v,d,Cd_inf,nu_f,rho_f);
        [Fdrag,phi,phi_flux] = calc_Fd_feedback_on_fluid(z,f_fluid,r,d,V,v,Lx,rho_f,u_star);
        find_bed_surface_i=1;
        while phi(find_bed_surface_i)>0.4
            find_bed_surface_i=find_bed_surface_i+1;
        end
        u=[zeros(find_bed_surface_i-1,1);u(1:end-find_bed_surface_i+1)];
        u_steady = [zeros(find_bed_surface_i-1,1);u_steady(1:end-find_bed_surface_i+1)];
        dz = z(2)-z(1);
        subplot(1,Nsubplots,Nflow)
        plot(u,z);
        hold on;
        plot(u/PARAMS.u_star,z,'b--');
        plot(u_steady,z,'k');
        plot(l/PARAMS.d0,zph,'r--');
        hold off;
        legend('u','u^+','l^+');
        title(['u(z), Shields = ',num2str(Shields),', u^* = ',num2str(u_star)]);
        xlabel('u [m/s]');
        ylabel('z [m]');
        ylim([0 Lz_plot]);
        xlim([-2 U_plot_scale]);
    end
    
    if FLAGS.saveresults > 0
        E         = zeros(It,1);
        E_height  = zeros(It,1);
        E_spring  = zeros(It,1);
        E_kin     = zeros(It,1);
        tvec      = zeros(It,1);
        Sediment_flux_time_series = zeros(It,1);
        phi_time_avg              = zeros(Iz-1,1);
        phi_flux_time_avg         = zeros(Iz-1,1);
        Fs = 1000; it_low_tres=1;
        r_vec_low_tres = zeros(3,P,round(It/round(1/Fs/dt))+1);
        v_vec_low_tres = zeros(3,P,round(It/round(1/Fs/dt))+1);
        t_vec_low_tres = zeros(round(It/round(1/Fs/dt))+1,1);
        u_low_tres     = zeros(Iz-1,round(It/round(1/Fs/dt))+1);
    end
    if FLAGS.saveresults > 1
        r_vec     = zeros(3,P,It);
        omega_vec = zeros(P,It);
        v_vec     = zeros(3,P,It);
    end
    
%     weight_bottom = 1;
%     weight_top    = 1;
%     if weight_bottom == weight_top
%         flow_weighting = ones(1,Iz-1)*weight_top;
%     else
%         flow_weighting = (weight_bottom:-(weight_bottom-weight_top)/(Iz-2):weight_top);
%     end
%     size(z)
%     size(flow_weighting)
%     pause;

    test_plot = 0;
    max_itter = 200;
    relax_vec   = zeros(max_itter,1);
    resnorm_vec = zeros(max_itter,1);
    
    time_counter = tic;
    for it=1:It
%         time_at_start(it)=toc(time_counter);
%         if it>1
%             disp(['ratio of time of flow calc = ',num2str( time_flow_calculation/(time_at_start(it)-time_at_start(it-1)) )]);
%         end

        %% calculate particle motions;
        list = update_list(list,r,d+2*a_crit,Lx);
        [delta,N,Int] = collision_detection(list,d,r,Lx,a_crit);
        [vc,vcn,vct,T] = calc_vc(v,d,omega,N,Int);
        delta_dot = calc_deltadot(vc,N);
        nu = calc_nu(vct);
        fn = forcen(delta,delta_dot,k,nu_n);
        ft = forcet(nu,nu_t,mu,fn);
        if V_bridge_star>0 % in the case of capilary bonds we have to check this:
            fn_cap = forcecap(delta,delta_dot,d,Int,V_bridge_star,theta_contact_angle);
            fn(delta<0)=-fn_cap(delta<0);
            ft(delta<0)=0;
        end
        %%
        if FLAGS.flow
                if it==1
                    for i_startup_flow=1:Istartup_flow
                        f_fluid = calc_fluid_forces(u,z,r,v,d,Cd_inf,nu_f,rho_f);
                        [Fdrag,phi,phi_flux] = calc_Fd_feedback_on_fluid(z,f_fluid,r,d,V,v,Lx,rho_f,u_star);
                        [u,l]=transient_lenght_scale_model_particle_interaction(u,l,dt,dz,nu_f,u_star,Fdrag/rho_f);
                    end
                end
                f_fluid = calc_fluid_forces(u,z,r,v,d,Cd_inf,nu_f,rho_f);
                [Fdrag,phi,phi_flux] = calc_Fd_feedback_on_fluid(z,f_fluid,r,d,V,v,Lx,rho_f,u_star);
%                 time_flow = tic;
                [u,l]=transient_lenght_scale_model_particle_interaction(u,l,dt,dz,nu_f,u_star,Fdrag/rho_f);
%                 time_flow_calculation = toc(time_flow);
                if FLAGS.saveresults > 0
                    phi_time_avg = phi_time_avg + phi/It;
                    phi_flux_time_avg = phi_flux_time_avg + phi_flux/It;
                    Sediment_flux_time_series(it) = sum(phi_flux)*dz;
                end
        else
            f_fluid = 0;
        end
        M  = calc_roll_resistance_momentum(Int,omega,mu_r,omega_max_r,vct);
    %     M(:)=0;

        %% calculate Energies of particles.
        [E(it),E_height(it),E_spring(it),E_kin(it)] = calc_E(r,v,m,delta,k,g);

        tvec(it)=t;
        if FLAGS.saveresults > 1
            omega_vec(:,it)=omega(3,:);
            r_vec(:,:,it)=r;  
            v_vec(:,:,it)=v;
        end
        
        
%         if it>1
%             if E(it)>E(it-1)*2
%                 test_plot=1;
%             end
%         end
        if FLAGS.saveresults > 0
            if mod(it-1,round(1/Fs/dt))==0
                r_vec_low_tres(:,:,it_low_tres) = r;
                v_vec_low_tres(:,:,it_low_tres) = v;
                t_vec_low_tres(it_low_tres)     = t;
                if FLAGS.flow
                    u_low_tres(:,it_low_tres)       = u;
                end
                it_low_tres=it_low_tres+1;
            end
        end
        
        if mod(it-1,round(It/(max(500*Tmax,50))))==0 && FLAGS.realtime_animation || test_plot
    %         figure(figpart);
            if FLAGS.showomega3 || FLAGS.flow
            subplot(1,Nsubplots,1);
            end
            hold off;
            for p=1:P
                myplot(r(:,p),d(p),'color',color(:,mod(p-1,7)+1));
                if p==1
                    hold on;
                end
                        if V_bridge_star>0
                        myplot(r(:,p),d(p)+2*a_crit,':','color',color(:,mod(p-1,7)+1));
                        end
                if FLAGS.plot_numbers; text(r(1,p),r(2,p),num2str(p)); end
                if r(1,p)-Lx+d(p)/2>0
                    myplot(r(:,p)-[Lx;0;0],d(p),'color',color(:,mod(p-1,7)+1));
                end
                if r(1,p)+Lx-d(p)/2<Lx
                    myplot(r(:,p)+[Lx;0;0],d(p),'color',color(:,mod(p-1,7)+1));
                end  
%                 myplot(r(:,p)+[Lx;0;0],d(p),'color',color(:,mod(p-1,7)+1));
%                 myplot(r(:,p)-[Lx;0;0],d(p),'color',color(:,mod(p-1,7)+1));
            end
            view(0,90);
            title(['time = ',num2str(t),' [s]']);
            xlabel('x');
            ylabel('y');
            zlabel('z');
    %         axis equal;
            xlim([0 Lx]);
            ylim([0 Lz_plot]);
            if FLAGS.showomega3
                subplot(1,Nsubplots,Nomega);
                for p=1:P
                    plot(tvec(it),omega(3,p),'color',color(:,mod(p-1,7)+1),'marker','.');
                    title(['v=',num2str(v(1,1)),' -omega*d/2=',num2str(-omega(3,1)*d(1)/2)]);
                    if p==1
                        hold on;
                    end
                end
                xlim([0 Tmax])
            end
            if FLAGS.flow
                subplot(1,Nsubplots,Nflow);
                plot(u,z);
                hold on;
%                 plot(l,zph);
                plot(v(1,:),r(2,:),'k.');
                plot(u/PARAMS.u_star,z,'b--');
                plot(l/PARAMS.d0,zph,'r--');
                plot(u_steady,z,'k');
                hold off;
                title(['u(z), Shields = ',num2str(Shields),', u^* = ',num2str(u_star)]);
                xlabel('u [m/s]');
                ylabel('z [m]');
                ylim([0 Lz_plot]);
                xlim([-2 U_plot_scale]);
                legend('u','u_p','u^+','l^+');
            end
            drawnow;
            if test_plot
                if it>1
                    r_dotdot
                    fn
                end
                pause; %#ok
            end
        end
    
        [r_dotdot,omega_dot] = calc_r_dotdot_omega_dot(fn,ft,f_fluid,m,d,I,N,T,ag,M,P,Int,Ldx,FLAGS);
        if it==1
            r_old = r;
            [r,v,omega] = timestep(r,v,omega,r_dotdot,omega_dot,dt,Lx,Ldx,FLAGS);
            r_current = r;
        else
            [r,v,omega] = timestep_verlet(r_old,r_current,omega,r_dotdot,omega_dot,dt,Lx,Ldx,FLAGS);
            r_old = r_current;
            r_current = r;
        end
        t=t+dt;
    end

    if FLAGS.plot_trajectories && FLAGS.saveresults > 1
        subplot(1,Nsubplots,1);
        for p=1:P
            plot3(squeeze(r_vec(1,p,:)),squeeze(r_vec(2,p,:)),squeeze(r_vec(3,p,:)),'color',color(:,mod(p-1,7)+1));
        end
    end

    if FLAGS.showEnergy
        figure;
        plot(tvec,E);
        hold on;
        plot(tvec,E_height);
        plot(tvec,E_spring);
        plot(tvec,E_kin);
        legend('E','E_{height}','E_{spring}','E_{kin}');
    end

    if FLAGS.showomega3_end_of_sim
        figure;
        plot(tvec,omega_vec);
        leg_vec = cell(P,1);
        for p=1:P; leg_vec{p}=num2str(p); end
        legend(leg_vec);
        title('angular velocity [rad/s]');
        xlabel('t');
        ylabel('omega [rad/s]');
    end

    %% save input data in a structure;
%     PARAMS.P=P;
%     PARAMS.d=d;
%     PARAMS.Lx=Lx;
%     PARAMS.Lz=Lz;
%     PARAMS.ywall = ywall;
%     PARAMS.g = g;
%     PARAMS.rho_f = rho_f;
%     PARAMS.nu_f  = nu_f;
%     PARAMS.Iz    = Iz;
%     PARAMS.u_star = u_star;
%     PARAMS.Cd_inf = Cd_inf;
%     PARAMS.rho_p = rho_p;
%     PARAMS.d0 = d0;
%     PARAMS.q = q;
%     PARAMS.k = k;
%     PARAMS.nu_n = nu_n;
%     PARAMS.nu_t = nu_t;
%     PARAMS.mu_r = mu_r;
%     PARAMS.omega_max_r =omega_max_r;
% 
%     NUMSET.dt=dt;
%     NUMSET.Tmax=Tmax;
%     NUMSET.It=It;

    if FLAGS.saveresults
        if exist(['RUNS/',store_name],'dir')~=7
            mkdir(['RUNS/',store_name]);
        end
        save(['RUNS/',store_name,'/PARAMS.mat'],'PARAMS');
        save(['RUNS/',store_name,'/NUMSET.mat'],'NUMSET');
        if FLAGS.saveresults > 1
            save(['RUNS/',store_name,'/r_vec.mat'],'r_vec');
            save(['RUNS/',store_name,'/omega_vec.mat'],'omega_vec');
            save(['RUNS/',store_name,'/v_vec.mat'],'v_vec');
        end
        if FLAGS.saveresults >0
            save(['RUNS/',store_name,'/phi_time_avg.mat'],'phi_time_avg');
            save(['RUNS/',store_name,'/phi_flux_time_avg.mat'],'phi_flux_time_avg');
            save(['RUNS/',store_name,'/Sediment_flux_time_series.mat'],'Sediment_flux_time_series');
            save(['RUNS/',store_name,'/r.mat'],'r');
            save(['RUNS/',store_name,'/v.mat'],'v');
            save(['RUNS/',store_name,'/omega.mat'],'omega');
            save(['RUNS/',store_name,'/r_vec_low_tres.mat'],'r_vec_low_tres');
            save(['RUNS/',store_name,'/v_vec_low_tres.mat'],'v_vec_low_tres');
            save(['RUNS/',store_name,'/t_vec_low_tres.mat'],'t_vec_low_tres');
            save(['RUNS/',store_name,'/u_low_tres.mat'],'u_low_tres');
        end
        save(['RUNS/',store_name,'/d.mat'],'d');
        save(['RUNS/',store_name,'/tvec.mat'],'tvec');
        save(['RUNS/',store_name,'/E.mat']       ,'E');
        save(['RUNS/',store_name,'/E_height.mat'],'E_height');
        save(['RUNS/',store_name,'/E_spring.mat'],'E_spring');
        save(['RUNS/',store_name,'/E_kin.mat']   ,'E_kin');
        if FLAGS.flow
            save(['RUNS/',store_name,'/z.mat']   ,'z');
            save(['RUNS/',store_name,'/u.mat']   ,'u');
        end
    end

end


function f_fluid = calc_fluid_forces(u,z,r,v,d,Cd_inf,nu_f,rho_f)
    P=numel(d);
    u_fluid_p =  zeros(3,P);
    dz = z(2)-z(1);
    z_p = r(2,:);
    for i=1:numel(z)-1
        for p=1:P
            if z(i+1)>z_p(p) && z(i)<z_p(p)
                u_fluid_p(1,p) = (u(i+1)-u(i))/dz*(z_p(p)-z(i)) + u(i);
            elseif z_p(p)>z(end)
                u_fluid_p(1,p) = v(1,p);
            end
        end
    end
    
    RuC = 24;
    Ru=zeros(P,1);
    f_fluid=zeros(3,P);
    for p=1:P
        norm_diff_u = norm(u_fluid_p(:,p)-v(:,p),2);
        Ru(p) = norm_diff_u*d(p)/nu_f;
        Cd = (sqrt(Cd_inf)+sqrt(RuC/Ru(p)))^2;
        if norm_diff_u==0
            Cd=1;
        end
        f_fluid(:,p) = pi/8*rho_f*d(p)^2*Cd*norm_diff_u*(u_fluid_p(:,p)-v(:,p));
    end
end

function [Fdrag,phi,phi_flux] = calc_Fd_feedback_on_fluid(z,f_fluid,r,d,V,v,Lx,rho_f,u_star)
    P=numel(d);
    Iz = numel(z);
    dz=z(2)-z(1);
    Fz = zeros(Iz,1);
    phi = zeros(Iz,1);
    phi_flux = zeros(Iz,1);
    tau_p_instant = zeros(Iz,1);  
    for p=1:P
        iz=1;
        R = d(p)/2;
        while r(2,p)-R>z(iz)+dz/2 && iz<Iz
            iz=iz+1;
        end
        if r(2,p)-R<z(iz)-dz/2
            x0 = z(iz)-r(2,p)-dz/2; % the case for particles fixed at the botom of the domain.
        else
            x0 = -R;
        end
        if r(2,p)+R>z(iz)+dz/2 && iz<Iz % added the part: && iz<Iz to prevent to get a large negative upperbound if the particle goes above the flow.
            x1 = z(iz)-r(2,p)+dz/2;
        else
            x1 = R;
        end
        A = Lx*d(p);
        Vpart   = pi*(R^2*x1-1/3*x1^3) - pi*(R^2*x0-1/3*x0^3);
        phi(iz)      = phi(iz)       + 1/(A*dz)*Vpart;                % duran et al. 2012 eq.(10)
        phi_flux(iz) = phi_flux(iz)  + 1/(A*dz)*Vpart*v(1,p); 
        Fz(iz)  = Fz(iz)  - 1/(A*dz)*f_fluid(1,p)*Vpart/V(p);   % duran et al. 2012 eq.(8)
        while r(2,p)+R>z(iz)+dz/2 && iz<Iz
            iz = iz+1;
            x0 = x1;
            if r(2,p)+R>z(iz)+dz/2
                x1 = z(iz)+dz/2 - r(2,p);
            else
                x1 = R;
            end
            Vpart   = pi*(R^2*x1-1/3*x1^3) - pi*(R^2*x0-1/3*x0^3);
            phi(iz)      = phi(iz)       + 1/(A*dz)*Vpart;
            phi_flux(iz) = phi_flux(iz)  + 1/(A*dz)*Vpart*v(1,p); 
            Fz(iz)  = Fz(iz)  - 1/(A*dz)*f_fluid(1,p)*Vpart/V(p);
        end
        if iz==Iz
            disp(['warning: particle ',num2str(p),' went above the modeled flow']);
        end
    end
%     Fdrag = Fz./phi;
    Fdrag = Fz./(1-phi);
%     Fdrag(Fz==0)=0;
end

function [tau_p_instant,phi] = calc_tau_p_instant(z,f_fluid,r,d,V,it,Lx,rho_f,u_star)
    P=numel(d);
    Iz = numel(z);
    dz=z(2)-z(1);
    Fz = zeros(Iz,1);
    phi = zeros(Iz,1);
    tau_p_instant = zeros(Iz,1);
    for p=1:P
        iz=1;
        R = d(p)/2;
        while r(2,p)-R>z(iz)+dz/2 && iz<Iz
            iz=iz+1;
        end
        if r(2,p)-R<z(iz)-dz/2
            x0 = z(iz)-r(2,p)-dz/2;
        else
            x0 = -R;
        end
        if r(2,p)+R>z(iz)+dz/2
            x1 = z(iz)-r(2,p)+dz/2;
        else
            x1 = R;
        end
        A = Lx*d(p);
        Vpart   = pi*(R^2*x1-1/3*x1^3) - pi*(R^2*x0-1/3*x0^3);
        phi(iz)      = phi(iz) + 1/(A*dz)*Vpart;                % duran et al. 2012 eq.(10)
        Fz(iz)  = Fz(iz)  - 1/(A*dz)*f_fluid(1,p)*Vpart/V(p);   % duran et al. 2012 eq.(8)
        while r(2,p)+R>z(iz)+dz/2 && iz<Iz
            iz = iz+1;
            x0 = x1;
            if r(2,p)+R>z(iz)+dz/2
                x1 = z(iz)+dz/2 - r(2,p);
            else
                x1 = R;
            end
            Vpart   = pi*(R^2*x1-1/3*x1^3) - pi*(R^2*x0-1/3*x0^3);
            phi(iz)      = phi(iz) + 1/(A*dz)*Vpart;
            Fz(iz)  = Fz(iz)  - 1/(A*dz)*f_fluid(1,p)*Vpart/V(p);
        end
        if iz==Iz
            disp(['warning: particle ',num2str(p),' went above the modeled flow']);
        end
    end
    Fdrag = Fz./phi;
    Fdrag(Fz==0)=0;
    tau_p_instant(Iz) = phi(Iz)/(1-phi(Iz))*Fdrag(Iz);
    for iz=Iz-1:-1:1
        tau_p_instant(iz) = tau_p_instant(iz+1) + phi(iz)/(1-phi(iz))*Fdrag(iz)*dz;  % duran et al. (2012), eq.(18)
%         tau_p_instant(iz) = tau_p_instant(iz+1) + 0.5*phi(iz)/(1-phi(iz))*Fdrag(iz)*dz + 0.5*phi(iz+1)/(1-phi(iz+1))*Fdrag(iz+1)*dz;  % duran et al. (2012), eq.(18)
    end
    tau_f = rho_f*u_star^2;
    tau_p_instant = - tau_p_instant;
    
    
%     if tau_p_instant(1)>rho_f*u_star^2
%         tau_p_instant = tau_p_instant/tau_p_instant(1)*rho_f*u_star^2;
%     end
    
%     figure;
%     plot(tau_p_instant,z);
%     hold on;
%     plot(tau_f*ones(size(z)),z);
%     pause;
%     tau_p_instant = min(tau_p_instant,tau_f);
end

function [tau_p_instant,phi] = calc_tau_p_instant_BACKUP(z,f_fluid,r,d,V,it,A,rho_f,u_star)
    P=numel(d);
    Iz = numel(z);
    dz=z(2)-z(1);
    Fz = zeros(Iz,1);
    phi = zeros(Iz,1);
    tau_p_instant = zeros(Iz,1);
    for p=1:P
        %%
        %{
        iz=1;
        R = d(p)/2;
        while r(2,p)-R>z(iz) && iz<Iz
            iz=iz+1;
        end
        istart = 1;
        while r(2,p)+R>z(iz) && iz<Iz
            if iz>1 && iz<=Iz
                if istart == 1
                    if z(iz-1)==0
                        x0=0;
                    else
                        x0 = -R;
                    end
                    x1 = z(iz)-r(2,p);
                else
                    x0 = z(iz-1)-r(2,p);
                    x1 = z(iz)  -r(2,p);
                end
                istart=istart+1;
                Vpart   = pi*(R^2*x1-1/3*x1^3) - pi*(R^2*x0-1/3*x0^3);
                phi(iz) = phi(iz) + 1/(A*dz)*Vpart;
                Fz(iz)  = Fz(iz)  - 1/(A*dz)*f_fluid(1,p)*Vpart/V(p);
            end
            iz=iz+1;
        end
        x0 = z(iz)-dz-r(2,p);
        x1 = R;
        Vpart   = pi*(R^2*x1-1/3*x1^3) - pi*(R^2*x0-1/3*x0^3);
        phi(iz) = phi(iz) + 1/(A*dz)*Vpart;
        Fz(iz)  = Fz(iz)  - 1/(A*dz)*f_fluid(1,p)*Vpart/V(p);
        if iz==Iz
            disp(['warning: particle ',num2str(p),' went above the modeled flow']);
        end
        %}
        %%
        iz=1;
        R = d(p)/2;
        while r(2,p)-R>z(iz)+dz/2 && iz<Iz
            iz=iz+1;
        end
        if r(2,p)-R<z(iz)-dz/2
            x0 = z(iz)-r(2,p)-dz/2;
        else
            x0 = -R;
        end
        if r(2,p)+R>z(iz)+dz/2
            x1 = z(iz)-r(2,p)+dz/2;
        else
            x1 = R;
        end
        Vpart   = pi*(R^2*x1-1/3*x1^3) - pi*(R^2*x0-1/3*x0^3);
        phi(iz)      = phi(iz) + 1/(A*dz)*Vpart;
        Fz(iz)  = Fz(iz)  - 1/(A*dz)*f_fluid(1,p)*Vpart/V(p);
        while r(2,p)+R>z(iz)+dz/2 && iz<Iz
            iz = iz+1;
            x0 = x1;
            if r(2,p)+R>z(iz)+dz/2
                x1 = z(iz)+dz/2 - r(2,p);
            else
                x1 = R;
            end
            Vpart   = pi*(R^2*x1-1/3*x1^3) - pi*(R^2*x0-1/3*x0^3);
            phi(iz)      = phi(iz) + 1/(A*dz)*Vpart;
            Fz(iz)  = Fz(iz)  - 1/(A*dz)*f_fluid(1,p)*Vpart/V(p);
        end
        if iz==Iz
            disp(['warning: particle ',num2str(p),' went above the modeled flow']);
        end
    end
    Fdrag = Fz./phi;
    Fdrag(Fz==0)=0;
    tau_p_instant(Iz) = phi(Iz)/(1-phi(Iz))*Fdrag(Iz);
    for iz=Iz-1:-1:1
        tau_p_instant(iz) = tau_p_instant(iz+1) + phi(iz)/(1-phi(iz))*Fdrag(iz)*dz;
    end
    tau_f = rho_f*u_star^2;
    tau_p_instant = - tau_p_instant;
    
    
%     if tau_p_instant(1)>rho_f*u_star^2
%         tau_p_instant = tau_p_instant/tau_p_instant(1)*rho_f*u_star^2;
%     end
    
%     figure;
%     plot(tau_p_instant,z);
%     hold on;
%     plot(tau_f*ones(size(z)),z);
%     pause;
%     tau_p_instant = min(tau_p_instant,tau_f);
end


function phi=redistribute_phi(phi,phi_max)
    phi_excess = 0;
    for i=1:numel(phi)
        if phi(i)>phi_max
            phi_diff = phi(i)-phi_max;
            phi_excess = phi_excess + phi_diff;
            phi(i) = phi_max;
        elseif phi(i)<phi_max
            phi_diff = phi_max - phi(i);
            if phi_diff<phi_excess
                phi(i)=phi_max;
                phi_excess = phi_excess - phi_diff;
            else
                phi(i)=phi(i)+phi_excess;
                phi_excess = 0;
            end
        end
    end
end

function [tau_p,tau_p_instant] = calc_tau_p(z,tau_p,f_fluid,r,d,V,it,A,rho_f,u_star)
    P=numel(d);
    Iz = numel(z);
    dz=z(2)-z(1);
    Fz = zeros(Iz,1);
    phi = zeros(Iz,1);
    tau_p_instant = zeros(Iz,1);
    for p=1:P
        iz=1;
        while r(2,p)>=z(iz) && iz<Iz
            iz=iz+1;
        end
        Fz(iz) = Fz(iz) - 1/(A*dz)*f_fluid(1,p);
        phi(iz) = phi(iz) + 1/(A*dz)*V(p);
%         tau_p_instant(1:iz)=tau_p_instant(1:iz)+f_fluid(1,p)/(A*dz);
        if iz==Iz
            disp('warning: particle went above the modeled flow');
        end
    end
    Fdrag = Fz./phi;
    Fdrag(Fz==0)=0;
    tau_p_instant(Iz) = phi(Iz)/(1-phi(Iz))*Fdrag(Iz);
    for iz=Iz-1:-1:1
        tau_p_instant(iz) = tau_p_instant(iz+1) + phi(iz)/(1-phi(iz))*Fdrag(iz)*dz;
    end
%     figure;
%     plot(tau_p_instant);
%     pause;
    tau_f = rho_f*u_star^2;
    tau_p_instant=min(tau_p_instant,tau_f);
    
%     figure;
%     plot(tau_p);
%     tau_f
%     pause;
    disp('hier ben ik gebleven');
    pause;
%     tau_p_instant=tau_p_instant/max(tau_p_instant)*tau_f;
    tau_p = 1/it*( tau_p_instant + (it-1)*tau_p);
%     if it<=3
%         tau_p = tau_p_instant;
%     end
        
end

function [r,v,omega]=timestep(r,v,omega,r_dotdot,omega_dot,dt,Lx,Ldx,FLAGS)
    [~,I]=size(r_dotdot);
    for i=1:I
        if i>Ldx || FLAGS.particles_fixed_at_bed==0
            omega(:,i) = omega(:,i) + dt*omega_dot(:,i);
            r(:,i) = r(:,i)+dt*v(:,i);
            r(1,i) = mod(r(1,i),Lx);
            v(:,i) = v(:,i)+dt*r_dotdot(:,i);
        end
    end
end

function [r,v,omega] = timestep_verlet(r_old,r_current,omega,r_dotdot,omega_dot,dt,Lx,Ldx,FLAGS)
    r=r_old;
%     if FLAGS.particles_fixed_at_bed
%         r(:,Ldx+1:end) = 2*r_current(:,Ldx+1:end) - r_old(:,Ldx+1:end) + dt^2*r_dotdot(:,Ldx+1:end);
%     else
        r = 2*r_current - r_old + dt^2*r_dotdot;
%     end
    r(1,:)=mod(r(1,:),Lx);
    v_old = compute_velocity(r,r_old,dt,Lx); %(r-r_old)/(2*dt);
    v     = v_old + dt*r_dotdot;
    omega = omega + dt*omega_dot;
    if FLAGS.particles_fixed_at_bed
        v(:,1:Ldx)=0;
        omega(:,1:Ldx)=0;
    end
end

function v_old = compute_velocity(r,r_old,dt,Lx)
    index = abs(r(1,:)-r_old(1,:))>Lx/4;
    r_periodic = r(1,index);
    r_periodic_old = r_old(1,index);
    ind = r_periodic-r_periodic_old<0;
    r_periodic(ind)=r_periodic(ind)+Lx;
    r_periodic(ind~=1)=r_periodic(ind~=1)-Lx;
    r(1,index)=r_periodic;
    v_old = (r-r_old)/(2*dt);
end

function [E,E_height,E_spring,E_kin] = calc_E(r,v,m,delta,k,g)
    E_height = sum(r(2,:).*m*g);
    E_spring = sum(0.5*delta.^2*k);
    E_kin = sum(0.5*m.*sum(v.*v,1));
    E = E_height + E_spring + E_kin;
end

function [r_dotdot,omega_dot] = calc_r_dotdot_omega_dot(fn,ft,f_fluid,m,d,I,N,T,ag,M,P,Int,Idx,FLAGS)
    r_dotdot  = repmat(ag,1,P);
    omega_dot = zeros(3,P);
    [~,int]=size(Int);
    for i=1:int
        r_dotdot(:,Int(1,i))      = r_dotdot(:,Int(1,i))  - (fn(i)*N(:,i) + ft(i)*T(:,i))/m(Int(1,i));
        omega_dot(:,Int(1,i))     = omega_dot(:,Int(1,i)) + cross(-d(Int(1,i))/2*N(:,i),ft(i)*T(:,i))/I(Int(1,i)) + M(:,i)/I(Int(1,i));
        if Int(2,i)~=0
            r_dotdot(:,Int(2,i))  = r_dotdot(:,Int(2,i))  + (fn(i)*N(:,i) + ft(i)*T(:,i))/m(Int(2,i));
            omega_dot(:,Int(2,i)) = omega_dot(:,Int(2,i)) + cross(-d(Int(2,i))/2*N(:,i),ft(i)*T(:,i))/I(Int(2,i)) + M(:,i+int)/I(Int(2,i));
        end
    end
    r_dotdot = r_dotdot+f_fluid./m;
    if FLAGS.particles_fixed_at_bed 
        r_dotdot(:,1:Idx)=0;
        omega_dot(:,1:Idx)=0;
    end
end

function [fn]=forcen(delta,deltadot,k,nu_n)
    fn = k*delta+nu_n*deltadot;
end

function [ft]=forcet(nu,nu_t,mu,fn)
    %% old
%     ft = -nu_t*nu;
    %% new; 
    ft = -sign(nu).*min(abs(nu_t*nu),mu*fn); % viscouse friction proposed in equation (48) in Luding (1998)
end

function f_cap = forcecap(delta,delta_dot,d,Int,V_bridge_star,theta)
%     Int
%     size(Int)
    I=numel(delta);
    R = zeros(size(delta));
    for i=1:I
        if Int(2,i)==0
            R(i) = d(Int(1,i))/2;
        else
            R(i) = d(Int(1,i)).*d(Int(2,i))./(d(Int(1,i))+d(Int(2,i)));
        end
    end
%     R = d(Int(1,:)).*d(Int(2,:))./(d(Int(1,:))+d(Int(2,:)));
    gamma = 20.6e-3; % [N/m]; table 1 in Gladkyy & Schwarze (2014): 20.6e-3 [N/m];
    a = -delta;
    Ca     = 1+6*a./(2*R); 
    Ctheta = 1+1.1*sin(theta); 
    V_bridge = V_bridge_star.*R.^3; 
    beta = asin( (V_bridge./(0.12*(2*R).^3.*Ca*Ctheta)).^(0.25)  );
%     a_crit_star = (1+theta/2)*(V_bridge_star^(1/3)+V_bridge_star^(2/3)/10);
%     a_crit = a_crit_star*R
%     size(R)
    R1 = (R.*(1-cos(beta))+a)./cos(beta+theta);
    R2 = R.*sin(beta)+R1.*(sin(beta+theta)-1);
    pk = gamma*(1./R1+1./R2);
    f_cap = pi/4*(2*R).^2.*sin(beta).^2.*pk + gamma*pi*2*R.*sin(beta).*sin(beta+theta);
    f_cap(delta_dot>0)=0;
end

function M  = calc_roll_resistance_momentum(Int,omega,mu_r,omega_max_r,vct)
    [~,I]=size(Int);
    M=zeros(3,2*I);
    for i=1:I
        M(:,i)   = - sign(omega(:,Int(1,i))).*min(abs(omega(:,Int(1,i))),omega_max_r)*mu_r;
        if Int(2,i)~=0
            M(:,i+I) = - sign(omega(:,Int(2,i))).*min(abs(omega(:,Int(2,i))),omega_max_r)*mu_r;
        end
    end
end

function nu = calc_nu(vct)
    nu = sqrt(sum(vct.*vct,1));
end

function deltadot = calc_deltadot(vc,N)
    deltadot = sum(-vc.*N,1);
end

function [vc,vcn,vct,T] = calc_vc(v,d,omega,N,Int)
    [~,I]=size(Int);
    vc  = zeros(3,I);
    vcn = zeros(3,I);
    vct = zeros(3,I);
    T   = zeros(3,I);
    for i=1:I
        if Int(2,i)==0
            vc (:,i) =              -v(:,Int(1,i))-cross(d(Int(1,i))/2*omega(:,Int(1,i)),N(:,i));
        else
            vc (:,i) = v(:,Int(2,i))-v(:,Int(1,i))-cross(d(Int(1,i))/2*omega(:,Int(1,i))+d(Int(2,i))/2*omega(:,Int(2,i)),N(:,i));
        end
        vcn(:,i) = N(:,i)*(vc(:,i)'*N(:,i));
        vct(:,i) = vc(:,i)-vcn(:,i);
        if norm(vct(:,i),2)==0
            T(:,i)  = rotation_matrix(pi/2)*N(:,i);
        else
            T(:,i)  = vct(:,i)/norm(vct(:,i),2);
        end
    end
end


function [be_list] = find_first_begining(be_list,P)
    i = 1;
    while be_list(i)>P
        i=i+1;
    end
    be_list = [ be_list(i:end); be_list(1:i-1)];
end

function [delta,N,Int] = collision_detection(list,d,r,Lx,a_crit)
    [I,~]=size(list);
    P = I/2;
    be_list = list(list(:,4),1);
    be_list = find_first_begining(be_list,P);
    ip =1;
    Ir=[];
    Ic=[];
    while ip<=I
        [Ir,Ic,~] = find_neighbours_xy(be_list,list,Ir,Ic,ip,P);
        ip=ip+1;
        if ip>I
            break
        end
        while be_list(ip)>P % skip end points of particles;
            ip=ip+1;
            if ip>I
                break
            end
        end
    end
    % up till now these are colliding bounding boxes; particles dont have to collide.    
    %% interparticle interactions;
    delta = zeros(1,numel(Ir));
    N     = zeros(3,numel(Ir));
    for i=1:numel(Ir)
        [delta(i),N(:,i)] = calc_delta(d(Ir(i)),d(Ic(i)),r(:,Ir(i)),r(:,Ic(i)),Lx);
    end
    positive_delta = delta>-a_crit;  % without moisture effects a_crit = 0;
    N = N(:,positive_delta);
    delta = delta(positive_delta);
    Int = [Ic(positive_delta);Ir(positive_delta)];
    %% wall interactions;
    
    positive_delta_wall = list(r(2,:)<d/2,1)';
    N     = [N    , repmat([0;-1; 0],1,numel(positive_delta_wall))]; 
    delta = [delta, d(positive_delta_wall)/2-r(2,positive_delta_wall)];
    Int   = [Int, [positive_delta_wall;zeros(1,numel(positive_delta_wall))]];
end
    
function [Ir,Ic,ip_next] = find_neighbours_x(be_list,list,Ir,Ic,ip,P)
    ipc=mod(ip,2*P)+1;
    first_other = 0;
    ip_next = ipc;
    while be_list(ipc)~=be_list(ip)+P
        if be_list(ipc)<=P
            Ir = [Ir be_list(ip) ];
            Ic = [Ic be_list(ipc)];
            if first_other==0
                first_other = 1;
                ip_next = ipc;
            end
        end
        ipc=mod(ipc,2*P)+1;
    end
end

function [Ir,Ic,ip_next] = find_neighbours_xy(be_list,list,Ir,Ic,ip,P)
    ipc=mod(ip,2*P)+1;
    first_other = 0;
    ip_next = ipc;
    while be_list(ipc)~=be_list(ip)+P
        if be_list(ipc)<=P  %% based on x-dir this is a possible colliding bounding box.
            %% now check the y-dir;
            be_list_y = [ ... 
                list(be_list(ip),3),...
                list(be_list(ipc),3),...
                list(be_list(ip)+P,3),...
                list(be_list(ipc)+P,3)];
            [~,ind]=sort(be_list_y);
            if ind(1)+2~=ind(2)
                Ir = [Ir be_list(ip) ];
                Ic = [Ic be_list(ipc)];
            end
            if first_other==0
                first_other = 1;
                ip_next = ipc;
            end
        end
        ipc=mod(ipc,2*P)+1;
    end
end

function plot_bounding_boxes(list,color,Lx,Lz)
    [I,~]=size(list);
    P=I/2;
    figure;
    ax(1,2) = subplot(2,2,2);
        hold on;
        for p=1:P
            if list(p,2)<list(p+P,2)
                plot([list(p,2) list(p+P,2) list(p+P,2) list(p,2) list(p,2)],[list(p,3) list(p,3) list(p+P,3) list(p+P,3) list(p,3)],'color',color(:,mod(p-1,7)+1));
                text( 0.5*(list(p,2)+list(p+P,2)),0.5*(list(p,3)+list(p+P,3)),num2str(p));
            else
                plot([list(p,2)-Lx list(p+P,2)    list(p+P,2)    list(p,2)-Lx list(p,2)-Lx],[list(p,3) list(p,3) list(p+P,3) list(p+P,3) list(p,3)],'color',color(:,mod(p-1,7)+1));
                plot([list(p,2)    list(p+P,2)+Lx list(p+P,2)+Lx list(p,2)    list(p,2)   ],[list(p,3) list(p,3) list(p+P,3) list(p+P,3) list(p,3)],'color',color(:,mod(p-1,7)+1));
                text( 0.5*(list(p,2)-Lx+list(p+P,2)),0.5*(list(p,3)+list(p+P,3)),num2str(p));
                text( 0.5*(list(p,2)+list(p+P,2)+Lx),0.5*(list(p,3)+list(p+P,3)),num2str(p));
            end
        end
        plot([0 0]   ,[0 1*Lz],'k');
        plot([1 1]*Lx,[0 1*Lz],'k');
        grid on;
        title('bounding boxes');
    ax(1,1) = subplot(2,2,1); % y bounding boxes;
            hold on;
            for p=1:P
                plot([0 0]+p,[list(p,3) list(p+P,3)],'color',color(:,mod(p-1,7)+1));
            end
            xlim([0 P+1]);
    ax(2,2) = subplot(2,2,4); % x bounding boxes;
        hold on;
        for p=1:P
            if list(p,2)<list(p+P,2)
                plot([list(p,2) list(p+P,2)],[0 0]+p,'color',color(:,mod(p-1,7)+1));
            else
                plot([list(p,2)-Lx list(p+P,2)],[0 0]+p,'color',color(:,mod(p-1,7)+1));
                plot([list(p,2) list(p+P,2)+Lx],[0 0]+p,'color',color(:,mod(p-1,7)+1));
            end
        end
        plot([0 0]   ,[-1 P+2],'k');
        plot([1 1]*Lx,[-1 P+2],'k');
        ylim([0 P+1]);
        
        xmargl          = 0.1;
        xsubmarg        = 0.1;
        xwidthsmall     = 0.1;
        xmargr          = 0.1;
        ymargb          = 0.1;
        ymargt          = 0.1;
        ysubmarg        = 0.1;
        ywidthsmall     = 0.1;
        xwidthlarge     = 1-xmargl-xmargr-xsubmarg-xwidthsmall;
        ywidthlarge     = 1-ymargb-ymargt-ysubmarg-ywidthsmall;
        for r = 1:2
            for c = 1:2
                xll = xmargl + (xwidthsmall+xsubmarg)*(c-1);
                yll = ymargb + (ywidthsmall+ysubmarg)*(2-r);
                if r==1 && c==1
                    xwidth = xwidthsmall;
                    ywidth = ywidthlarge;
                elseif r==1 && c==2
                    xwidth = xwidthlarge;
                    ywidth = ywidthlarge;
                else
                    xwidth = xwidthlarge;
                    ywidth = ywidthsmall;
                end
                if (r==2 && c==1)~=1
                    set(ax(r,c),'position',[xll yll xwidth ywidth]);
                end
            end
        end
        drawnow;
end

function list = initiate_list(r,d,Lx)
    [~,P] = size(r); 
    list = zeros(2*P,5);
    list(:,1)=1:2*P;
    list(1  :P  ,2) = r(1,:)-d/2; % x-dir;
    list(1  :P  ,3) = r(2,:)-d/2; % y-dir;
    list(P+1:2*P,2) = r(1,:)+d/2; % x-dir;
    list(P+1:2*P,3) = r(2,:)+d/2; % y-dir;
        list(list(:,2)<0 ,2) = list(list(:,2)<0 ,2) + Lx;
        list(list(:,2)>Lx,2) = list(list(:,2)>Lx,2) - Lx;
    [~,list(1:2*P,4)] = sort(list(:,2));
    [~,list(1:2*P,5)] = sort(list(:,3));
end

function list = update_list(list,r,d,Lx)
    [~,P] = size(r); 
    list(1  :P  ,2) = r(1,:)-d/2; % x-dir;
    list(1  :P  ,3) = r(2,:)-d/2; % y-dir;
    list(P+1:2*P,2) = r(1,:)+d/2; % x-dir;
    list(P+1:2*P,3) = r(2,:)+d/2; % y-dir;
        list(list(:,2)<0 ,2) = list(list(:,2)<0 ,2) + Lx;
        list(list(:,2)>Lx,2) = list(list(:,2)>Lx,2) - Lx;
    [~,indx] = sort(list(list(:,4),2));  % a very well sorted list already;
    [~,indy] = sort(list(list(:,5),3));  % a very well sorted list already;
    list(:,4)=list(indx,4);
    list(:,5)=list(indy,5);
end

% function [delta,N] = calc_delta(d1,d2,r1,r2,L)
%     if abs(r1(1)-r2(1))<L/2
%     elseif r1(1)-r2(1)<0
%         r2(1) = r2(1)-L;
%     elseif r1(1)-r2(1)>0
%         r2(1) = r2(1)+L;
%     else
%         disp('error, I did not anticipate for this! (calc_delta)');
%     end
%     N = (r1-r2)/norm(r1-r2,2);
%     dr =(r1-r2)'*N;
%     delta = 0.5*(d1+d2)-dr;
% end

function s = myplot3(x,d,varargin)
    [X,Y,Z]=sphere(10);
    s = surf(0.5*d*X+x(1),0.5*d*Y+x(2),0.5*d*Z+x(3),varargin{:});       
end

function s = myplot(x,d,varargin)
    theta = linspace(0,2*pi,20);
    X = d/2*cos(theta);
    Y = d/2*sin(theta);
    s = plot(x(1)+X,x(2)+Y,varargin{:});      
end

function rot = rotation_matrix(t)
rot = [ cos(t) -sin(t) 0; 
        sin(t)  cos(t) 0;
        0       0      1];
end

function [L2]=L2norm(u)
    I=length(u(:));
    L2=sqrt(sum(u(:).*conj(u(:)))/I);
end
