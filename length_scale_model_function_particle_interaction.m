function [u,z,l,zph,tau_f,tau_f_nu,tau_f_t]=length_scale_model_function_particle_interaction(rho_f,nu,L,Iz,u_star,tau_p,d0,phi)
    

    % parameters;
    % rho_f = 1.225;
    % nu = 1.5e-5;
    RvD = 26;
    Rc = 7;
    kappa = 0.41;
    phi_b = 0.6;
    lambda_b = d0*sqrt((1-phi_b)/(18*phi_b));
    
    tau_f = rho_f*u_star^2-tau_p;
%     alpha = 0.5;
%     tau_f = rho_f*u_star^2-(alpha*tau_p+(1-alpha)*[tau_p(2:end);0]);

    tau_f(tau_f<0)=0;
%     Th = 0.0001;
%     tau_f(tau_f<-Th*rho_f*u_star^2)=-Th*rho_f*u_star^2;
    
        
    z = linspace(0,L,Iz);
    if max(phi)<phi_b/2
        z_surface_bed = 0;
        ind_bed = 0;
    else
        logical=phi>phi_b/2;
        ind_bed = max(find(logical==1));
        z_surface_bed = z(ind_bed);
    end
    ind_bed_margin = 5;
    
    dz = z(2)-z(1);
    z=z(1:end-1);
    zph = z+0.5*dz;
    
    l = zeros(Iz-1,1);
    
%     figure;
%     plot(tau_f,z);
%     pause;
    
%     if dz>0.2*d
%         disp(['dz > 0.2*d  (',num2str(dz),' > ',num2str(0.2*d)]);
% %         pause;
%     end
%     if dz>nu/u_star
%         disp(['dz > nu/u_star  (',num2str(dz),' > ',num2str(nu/u_star)]);
% %         pause;
%     end
    
%     l(1) = kappa*dz/2*(1-exp(-1/RvD*dz/2*u_star/nu));
    
    
    
    max_itter=600; 
%     tol = 1e-10; 
    relax = 0.8;
    
%     tol_rel_u = 5e-6;
%     tol_abs_u = 1e-12;%1e-13;
%     tol_rel_l = 5e-6;
%     tol_abs_l = 1e-12;
    
    tol_rel_u = 5e-5;
    tol_abs_u = 1e-8;%1e-13;
    tol_rel_l = 5e-5;
    tol_abs_l = 1e-8;
    
    u=zeros(Iz-1,1);
    u(1)=eps;
    itter = 1;
    res = 1;
    tol_exponential_bed = 1e-14;
    
    if ind_bed ==0
        u(1) = 0;
        l(1) = kappa*dz/2*(1-exp(-1/RvD*dz/2*u_star/nu));
    else
        u(1) = 0;%-1e-3;%1e-3;
        l(1) = find_l1(1,dz,Rc,u_star,nu,u(1),rho_f,tau_f,kappa,RvD,lambda_b,tol_rel_u,tol_abs_u,tol_rel_l,tol_abs_l,max_itter,relax);
%         Iz_tmp = ind_bed;
%         [u,l] = integrate_flow(rho_f,nu,Rc,kappa,l,u,tau_f,Iz_tmp+10,dz,tol_rel_u,tol_abs_u,tol_rel_l,tol_abs_l,max_itter,relax);
%         res = u(1) - u(ind_bed-ind_bed_margin)*exp(-dz*(ind_bed-ind_bed_margin-1)/lambda_b);
%         figure;
%         semilogx(u(1:ind_bed+9),z(1:ind_bed+9));
%         hold on;
%         semilogx(l(1:ind_bed+9),z(1:ind_bed+9));
%         us = u(ind_bed-1);
%         semilogx(us*exp((z(1:ind_bed)-z_surface_bed)/lambda_b),z(1:ind_bed));
%         legend('u(z)','l(z)','u_s*exp(z/lambda_b');
%         l(1:20)
%         pause;
%         while abs(res)>tol_exponential_bed && itter<max_itter
%             l(1) = find_l1(dz,Rc,u_star,nu,u(1),rho_f,tau_f,kappa,RvD,lambda_b,tol_l,tol_u,max_itter,relax);
%             
%             ind_bed
%             
%             res_old = res;
%             %%
% 
%             res = u(1) - u(ind_bed-ind_bed_margin)*exp(-dz*(ind_bed-ind_bed_margin-1)/lambda_b);
%             if itter ==1
%                 dresdu1 = 1;
%             else
%                 dresdu1 = (res-res_old)/du1;
%             end
%             du1 = -relax*res/dresdu1;
%             u(1) = u(1) + du1;
%             itter = itter+1;
%             
%         end
%         if itter==max_itter
%             res
%             disp('exponential bed itteration did not converge! (length_scale_model, line 48)');
%             pause;
%         end
        
    end
    
    
    [u,l,tau_f_nu,tau_f_t] = integrate_flow(rho_f,nu,Rc,kappa,l,u,tau_f,Iz,dz,tol_rel_u,tol_abs_u,tol_rel_l,tol_abs_l,max_itter,relax);
    
%     figure;
%     subplot(2,2,1);
%         semilogx(u,z,'.-');
%         title('u(z)');
%     subplot(2,2,2);
%         semilogx(l,zph,'.-');
% %         hold on;
% %         semilogx(l(1:iz_start-1),zph(1:iz_start-1),'.-');
%         title('l(z)');
%     subplot(2,2,3);
%         plot(u,z,'.-');
% %         hold on;
% %         plot(u(1:iz_start-1),z(1:iz_start-1),'.-')
%         title('u(z)');
%     subplot(2,2,4);
%         plot(l,zph,'.-');
% %         hold on;
% %         plot(l(1:iz_start-1),zph(1:iz_start-1),'.-');
%         title('l(z)');

end

function [u,l,tau_f_nu,tau_f_t] = integrate_flow(rho_f,nu,Rc,kappa,l,u,tau_f,Iz,dz,tol_rel_u,tol_abs_u,tol_rel_l,tol_abs_l,max_itter,relax)
    tau_f_nu = zeros(size(tau_f));
    tau_f_t  = zeros(size(tau_f));
    for iz = 2:Iz-1
        %% u(iz) itteration:
        [u(iz),tau_f_nu(iz),tau_f_t(iz)] = itter_u(rho_f,nu,iz,l(iz-1),u(iz-1),tau_f(iz),dz,tol_rel_u,tol_abs_u,max_itter,relax);

        %% l(iz) itteration:
        l(iz) = itter_l(l(iz-1),u(iz),iz,nu,Rc,kappa,dz,tol_rel_l,tol_abs_l,max_itter,relax);
    end
end


function l = find_l1(iz,dz,Rc,u_star,nu,u1,rho_f,tau_f,kappa,RvD,lambda_b,tol_rel_u,tol_abs_u,tol_rel_l,tol_abs_l,max_itter,relax)
    l = kappa*dz/2*(1-exp(-1/RvD*dz/2*u_star/nu));
    u2 = itter_u(rho_f,nu,iz,l,u1,tau_f(1),dz,tol_rel_u,tol_abs_u,max_itter,relax);
    
%     tol = 1e-14;
    itter = 1;
%     res_abs = l - kappa^2*lambda_b^2/nu/Rc*(u1+u2)/2;
    res_abs = l - kappa^2*lambda_b^2/nu/Rc*abs(u1+u2)/2;
    while abs(res_abs)>tol_abs_l && itter<max_itter
        if itter ==1
            dresdl = 1;
        else
            dresdl = (res_abs - res_old)/dl;
        end    
        dl = - relax*res_abs/dresdl;
        l = l + dl;
        u2 = itter_u(rho_f,nu,iz,l,u1,tau_f(1),dz,tol_rel_u,tol_abs_u,max_itter,relax);
        res_old = res_abs;
%         res_abs = l - kappa^2*lambda_b^2/nu/Rc*(u1+u2)/2;
        res_abs = l - kappa^2*lambda_b^2/nu/Rc*abs(u1+u2)/2;
        itter = itter+1;
    end
end

function [u_i,tau_f_nu_i,tau_f_t_i] = itter_u(rho_f,nu,iz,l_im1,u_im1,tau_f_i,dz,tol_rel_u,tol_abs_u,max_itter,relax)
    u_i = u_im1+eps; % initial guess;
    itter = 0; res_abs_u = 1; res_rel_u = 1;
    tol_rel = 1e-6;
    % while abs(res)>tol && itter<max_itter
    while abs(res_rel_u)>tol_rel_u && abs(res_abs_u)>tol_abs_u && itter<max_itter
        itter=itter+1;
        res_abs_u = rho_f*(nu+l_im1^2*abs((u_i-u_im1)/dz))*(u_i-u_im1)/dz - tau_f_i;
        %% new;
        res_nu = rho_f*(nu)*(u_i-u_im1)/dz - tau_f_i;
        res_l2 = rho_f*(l_im1^2*abs((u_i-u_im1)/dz))*(u_i-u_im1)/dz - tau_f_i;
        res_rel_u = res_abs_u/min(abs([res_nu,res_l2]));
        %
        dresdui = 1/dz*rho_f*(nu+l_im1^2*abs((u_i-u_im1)/dz)) + rho_f*l_im1^2*sign((u_i-u_im1)/dz)/dz; 
        u_i = u_i -relax*res_abs_u/dresdui;
    end
    tau_f_nu_i = rho_f*(nu                              )*(u_i-u_im1)/dz;
    tau_f_t_i  = rho_f*(    l_im1^2*abs((u_i-u_im1)/dz) )*(u_i-u_im1)/dz;
    if itter==max_itter
        res_abs_u
        res_rel_u
        disp(['warning! u(iz=',num2str(iz),') itteration stopped at max_itter (length_scale_model, line 181)']);
        pause;
    end
end

function l_i = itter_l(l_im1,u_i,iz,nu,Rc,kappa,dz,tol_rel_l,tol_abs_l,max_itter,relax)
    l_i = l_im1; % initial guess;
%     max_itter
    itter = 0; res_abs_l = 1; res_rel_l = 1;
    %while abs(res)>tol && itter<max_itter
    while abs(res_rel_l)>tol_rel_l && abs(res_abs_l)>tol_abs_l && itter<max_itter
        itter=itter+1;
%         res_abs_l = (l_im1-l_i)/dz + kappa*(1-exp(-sqrt(1/Rc*(u_i*0.5*(l_i+l_im1))/nu)));
        res_abs_l = (l_im1-l_i)/dz + kappa*(1-exp(-sqrt(1/Rc*(abs(u_i)*0.5*(l_i+l_im1))/nu)));
        %% new
        res_liz   = l_i/dz;
        res_lizm1 = l_im1/dz;
        res_kappa = kappa;
%         res_kappa_exp = kappa*(-exp(-sqrt(1/Rc*(u_i*0.5*(l_i+l_im1))/nu)));
        res_kappa_exp = kappa*(-exp(-sqrt(1/Rc*(abs(u_i)*0.5*(l_i+l_im1))/nu)));
        res_rel_l = res_abs_l/min(abs([res_liz,res_lizm1,res_kappa,res_kappa_exp]));
        % 
%         dresdli = -1/dz - kappa*exp(-sqrt(1/Rc*(u_i*0.5*(l_i+l_im1))/nu))  *   0.5*(1/Rc* (u_i*0.5*(l_i+l_im1) )/nu)^(-0.5)   *   1/Rc*u_i*0.5/nu;
        dresdli = -1/dz - kappa*exp(-sqrt(1/Rc*(abs(u_i)*0.5*(l_i+l_im1))/nu))  *   0.5*(1/Rc* (abs(u_i)*0.5*(l_i+l_im1) )/nu)^(-0.5)   *   1/Rc*abs(u_i)*0.5/nu;
        l_i = l_i - relax*res_abs_l/dresdli;
    end
    if itter==max_itter
        l_im1
        l_i
        u_i
        res_abs_l
        res_rel_l
        disp(['warning! l(',num2str(iz),') itteration stopped at max_itter (length_scale_model, line 211)']);
        pause;
    end
end