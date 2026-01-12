function [u,l]=transient_lenght_scale_model_particle_interaction(u,l,dt,dz,nu_f,u_star,Fd)
    % possibly make a switch between implicit and direct method (direct is
    % faster, if the timestep allows it).
%      [u,l] = timestep_ul(u,l,dt,dz,nu_f,u_star,Fd);
    [u,l] = semi_implicit_timestep_ul(u,l,dt,dz,nu_f,u_star,Fd);
end

function [u,l] = timestep_ul(u,l,dt,dz,nu_f,u_star,Fd)
    tau = [(nu_f+l(1:end-1).^2.*abs((u(2:end)-u(1:end-1))/dz)).*(u(2:end)-u(1:end-1))/dz;u_star^2];
    u(2:end) = u(2:end) + dt*( (tau(2:end)-tau(1:end-1))/dz + Fd(2:end)); 
    l = solve_l(u,l,dz,nu_f,u_star);
end

function [u,l] = semi_implicit_timestep_ul(u,l,dt,dz,nu_f,u_star,Fd)
    Iz=numel(u);
    T = dt*(nu_f+l(1:end-1).^2.*abs((u(2:end)-u(1:end-1))/dz))/dz^2;
    % allocate:
    i=zeros((Iz-1)*3,1);
    j=zeros((Iz-1)*3,1);
    k=zeros((Iz-1)*3,1);
    B=zeros(Iz,1);
    % u(1) = 0; 
    i(1)=1;
    j(1)=1;
    k(1)=1; % or any other constant if that is better for the condition number;
    B(1)=0;
    % u(internal);
    % diagonal; u(iz);
    i(2:Iz-1)=2:Iz-1;
    j(2:Iz-1)=2:Iz-1;
    k(2:Iz-1)=1 + T(2:end)+T(1:end-1);
    % diagnal+1; u(iz+1);
    i(Iz:2*Iz-3) = 2:Iz-1;
    j(Iz:2*Iz-3) = 3:Iz;
    k(Iz:2*Iz-3) = - T(2:end);
    % diagnal-1; u(iz-1);
    i(2*Iz-2:3*Iz-5) = 2:Iz-1;
    j(2*Iz-2:3*Iz-5) = 1:Iz-2;
    k(2*Iz-2:3*Iz-5) = - T(1:end-1);
    B(2:Iz-1)=u(2:end-1) + dt*Fd(2:Iz-1); % eventueel kunnen hier nog wel particle drag force termen in komen.
    % u(Iz)
    i(3*Iz-4:3*Iz-3)=Iz;
    j(3*Iz-4:3*Iz-3)=[Iz-1,Iz];
    % u(Iz-1)
    k(3*Iz-4)=-T(end);
    % u(Iz);
    k(3*Iz-3)=1 + T(end);
    B(Iz) = u(end) + dt*u_star^2/dz + dt*Fd(end);
    A=sparse(i,j,k);
    u=A\B;
    l = solve_l(u,l,dz,nu_f,u_star);
end

function [l]=solve_l(u,l,dz,nu_f,u_star)
    Iz=numel(u);
    l(1)=itter_l_first(dz,nu_f,u_star);
    for iz = 2:Iz
        l(iz)=itter_l(l(iz),l(iz-1),u(iz),dz,nu_f);
    end
end

function [l_i]=itter_l_first(dz,nu_f,u_star)
    RvD = 26;
    kappa = 0.41;
    z=0.5*dz;
    l_i = kappa*z*(1-exp(-1/RvD*z*u_star/nu_f));
end

function [l_i]=itter_l(l_i,l_im1,u_i,dz,nu_f)
    Rc = 7;
    kappa = 0.4;
    tol = 1e-6;
    relax = 0.7;
    Res    = kappa*(1-exp(-sqrt(1/Rc*(abs(u_i)*0.5*(l_i+l_im1))/nu_f)))-(l_i-l_im1)/dz;
    itter = 0; max_itter = 100;
    while abs(Res)>tol && itter<max_itter
        itter=itter+1;
        dResdl = -1/dz - kappa*exp(-sqrt(1/Rc*(abs(u_i)*0.5*(l_i+l_im1))/nu_f))  *   0.5*(1/Rc* (abs(u_i)*0.5*(l_i+l_im1) )/nu_f)^(-0.5)   *   1/Rc*abs(u_i)*0.5/nu_f;
        Res    = kappa*(1-exp(-sqrt(1/Rc*(abs(u_i)*0.5*(l_i+l_im1))/nu_f)))-(l_i-l_im1)/dz;
        Res_log(itter)=Res;
        l_i=l_i-Res/dResdl*relax;
    end
    if itter==max_itter
        disp('warning itteration stopped due to max_itter');
        fig=figure;
        plot(Res_log);
        pause;
        close(fig); 
        itter
        Res
    end
end