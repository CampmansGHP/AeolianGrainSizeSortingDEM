clear;
close all;

FLAG_export_figs = 0;

fs = 10;

% load_name = ['GEN_STATIC_BED_P1500_d50_0.25_std_phi_',num2str(0.5)];

% load_name = '';
load_name = 'example_transport';


% load_name = ['BATCH_GENinput_static_bed_P15000_d250mu_bed1000p_Shield',num2str(5/100,'%1.2f')];
% load_name = ['BATCH_GENinput_static_bed_P1500_d250mu_bed100p_Shield',num2str(5/100,'%1.2f')];
% load_name = ['BATCH_GENinput_static_bed_P1500_d250mu_bed100p_Shield',num2str(6/100,'%1.2f'),'_modified_init_flow'];
show_animation = 0;
plot_tails = 1;

load('default_colors.mat','color');
load(['RUNS/',load_name,'/phi_time_avg.mat'],'phi_time_avg');
load(['RUNS/',load_name,'/phi_flux_time_avg.mat'],'phi_flux_time_avg');
load(['RUNS/',load_name,'/Sediment_flux_time_series.mat'],'Sediment_flux_time_series');
load(['RUNS/',load_name,'/d.mat'],'d');
load(['RUNS/',load_name,'/E.mat'],'E');
load(['RUNS/',load_name,'/PARAMS.mat'],'PARAMS');
load(['RUNS/',load_name,'/NUMSET.mat'],'NUMSET');
load(['RUNS/',load_name,'/tvec.mat'],'tvec');
load(['RUNS/',load_name,'/z.mat'],'z');
load(['RUNS/',load_name,'/Sediment_flux_time_series.mat'],'Sediment_flux_time_series');

load(['RUNS/',load_name,'/t_vec_low_tres.mat'],'t_vec_low_tres');
load(['RUNS/',load_name,'/v_vec_low_tres.mat'],'v_vec_low_tres');
load(['RUNS/',load_name,'/r_vec_low_tres.mat'],'r_vec_low_tres');
load(['RUNS/',load_name,'/u_low_tres.mat'],'u_low_tres');

if sum(r_vec_low_tres(1,:,end))==0
    t_vec_low_tres = t_vec_low_tres(1:end-1);
    v_vec_low_tres = v_vec_low_tres(:,:,1:end-1);
    r_vec_low_tres = r_vec_low_tres(:,:,1:end-1);
    u_low_tres     = u_low_tres(:,1:end-1);
end

P=PARAMS.P;
Lx=PARAMS.Ldx*PARAMS.d0;
Itlowres = numel(t_vec_low_tres);
maxu = max(u_low_tres(:));


% Cplus = 4.5;
% figure;
% plot(u_low_tres(:,1),z);
% hold on;
% plot(u_loglaw(z,PARAMS,Cplus),z+3e-3);
% title('u(z) at t=0');

figure('units','centimeters','position',[2 2 8.5 8.5],'color',[1 1 1],'name','maximum particle height vs time');
    plot(t_vec_low_tres,max(squeeze(r_vec_low_tres(2,:,:))));
    title('max particle height in time');

figure('units','centimeters','position',[2 2 8.5 8.5],'color',[1 1 1],'name','maximum velocity of flow and particles');
    plot(t_vec_low_tres,max(squeeze((v_vec_low_tres(1,:,:).^2+v_vec_low_tres(2,:,:).^2).^0.5)));
    hold on;
    plot([0 5],[1 1]*PARAMS.d0/NUMSET.dt);
    plot(t_vec_low_tres,max(u_low_tres));
    legend('u_{p,max}','u_{stab,max}','u_{flow,max}');
    title('max particle velocity in time');
    xlabel('t [s]');
    ylabel('u_{p,max} [m/s]');

mean_d = mean(d(101:PARAMS.P));
std_d = std(d(101:PARAMS.P));
Nbins_d=30;
Nbins_z = 200;
z_max_bin = 0.06; %PARAMS.Lz/2;
start_t = 2.5;
[N,Xd]=hist(d(101:PARAMS.P),Nbins_d);

binbounds_z = linspace(0,z_max_bin,Nbins_z+1);
Xz = binbounds_z(1:end-1)+(binbounds_z(2)-binbounds_z(1))/2; % bincenters;

% if strcmp(load_name,'BATCH_GENinput_static_bed_P15000_d250mu_bed1000p_Shield0.05')==1
%     disp('test');
%     load(['RUNS/',load_name,'/Qperbin_struct.mat'],'Qperbin_struct');
% %     Xz = Qperbin_struct.Xz;
% %     Xd = Qperbin_struct.Xd;
%     Qperbin_dz = Qperbin_struct.Qperbin_dz;
%     %     Qperbin_struct=struct();
%     %     Qperbin_struct.Qperbin_dz=Qperbin_dz;
%     %     Qperbin_struct.Xz=Xz;
%     %     Qperbin_struct.Xd=Xd;
%     %     save(['RUNS/',load_name,'/Qperbin_struct.mat'],'Qperbin_struct');
%     
% else
%     
%    
    if exist(['RUNS/',load_name,'/Qperbin_struct.mat'],'file')==2
        load(['RUNS/',load_name,'/Qperbin_struct.mat'],'Qperbin_struct');
        if min(Qperbin_struct.Xz==Xz) && min(Qperbin_struct.Xd==Xd)
            Qperbin_dz = Qperbin_struct.Qperbin_dz;
        else
            disp('The bin bounds have changed between the saved data and the ones provided in this script.');
            disp('Would you like to recalculate the distribution within the new bin bounds?');
            disp('<< press enter to continue>> or << Ctrl+C >> to cancel');
            pause;
            clear Qperbin_struct
            Qperbin_dz = calc_Q_perbin_dz(d,PARAMS.Ldx,Nbins_d,Nbins_z,Xd,Xz,v_vec_low_tres,r_vec_low_tres,t_vec_low_tres,start_t,Lx);
            Qperbin_struct.Xz=Xz;
            Qperbin_struct.Xd=Xd;
            Qperbin_struct.Qperbin_dz=Qperbin_dz;
            save(['RUNS/',load_name,'/Qperbin_struct.mat'],'Qperbin_struct');
            disp('The binned results have been saved.');
        end
    else
        Qperbin_dz = calc_Q_perbin_dz(d,PARAMS.Ldx,Nbins_d,Nbins_z,Xd,Xz,v_vec_low_tres,r_vec_low_tres,t_vec_low_tres,start_t,Lx);
        Qperbin_struct.Xz=Xz;
        Qperbin_struct.Xd=Xd;
        Qperbin_struct.Qperbin_dz=Qperbin_dz;
        save(['RUNS/',load_name,'/Qperbin_struct.mat'],'Qperbin_struct');
        disp('The binned results have been saved.');
    end
    
Qperbin_d = calc_Q_perbin(d,PARAMS.Ldx,Nbins_d,Xd,v_vec_low_tres,t_vec_low_tres,start_t,Lx);

figure('units','centimeters','position',[2 2 30 10],'color',[1 1 1],'name','grain size distribution and flux per bin(d,z)');
subplot(1,4,1);
    plot(1:PARAMS.P-PARAMS.Ldx,d(PARAMS.Ldx+1:PARAMS.P),'k.');
    xlim([0 PARAMS.P-PARAMS.Ldx]);
    ylabel('grain diameter, d [m]');
    ylim([Xd(1)-(Xd(2)-Xd(1))/2 Xd(end)+(Xd(2)-Xd(1))/2]);
    xlabel('particle number [-]');
subplot(1,4,2);
%     plot(N/1400*Nbins_d/3,X);
    plot_binned_data(N/1400*Nbins_d/3,Xd);
    hold on;
    plot([0 0.1],[1 1]*mean_d,'k'); text(0.075,mean_d+0.17*std_d,'\mu');
    plot([0 0.1],[1 1]*mean_d+std_d,'k'); text(0.075,mean_d+(1+0.17)*std_d,'\mu+\sigma');
    plot([0 0.1],[1 1]*mean_d-std_d,'k'); text(0.075,mean_d+(-1+0.17)*std_d,'\mu-\sigma');
    xlabel('distribution [-]');
    ylabel('grain diameter, d [m]');
    ylim([Xd(1)-(Xd(2)-Xd(1))/2 Xd(end)+(Xd(2)-Xd(1))/2]);
    text(0.01,mean_d+2.5*std_d,['mean, \mu = ',num2str(mean_d,'%1.2e'),' [m]']);
    text(0.01,mean_d+2.1*std_d,['std, \sigma = ',num2str(std_d,'%1.2e'),' [m]']);
%     title(['mean, \mu=',num2str(mean_d,'%1.2e'),'[m]; std, \sigma = ',num2str(std_d,'%1.2e'),'[m]']);
subplot(1,4,3);
    
%     plot(Qperbin*PARAMS.rho_p*Nbins_d,X);
    plot_binned_data(Qperbin_d*PARAMS.rho_p*Nbins_d,Xd);
    xlabel('flux per bin [kg/m/s]');
    xl = xlim;
    xl(1)=0;
    xlim(xl);
    ylim([Xd(1)-(Xd(2)-Xd(1))/2 Xd(end)+(Xd(2)-Xd(1))/2]);
    title({'Q as function of grainsize bin';['Qtot = ',num2str(sum(Qperbin_d)*PARAMS.rho_p),' [kg/m/s]'];' '});  
subplot(1,4,4);
    [Xdz,Ydz]=ndgrid(Xd,Xz);
%     imagesc(Xd,Xz,Qperbin_dz*PARAMS.rho_p*Nbins_d*Nbins_z);
    my_tile_pcolor(Xz,Xd,Qperbin_dz'*PARAMS.rho_p*Nbins_d*Nbins_z);
%     pcolor(Ydz,Xdz,Qperbin_dz*PARAMS.rho_p*Nbins_d*Nbins_z);
    shading flat;
    xlim([0 z_max_bin]);
    ylim([Xd(1)-(Xd(2)-Xd(1))/2 Xd(end)+(Xd(2)-Xd(1))/2]);
    colorbar;
    colormap jet;
    xlabel('z [m]');
    title('flux per grainsize per height');
    

    
particle_velocity_sum_all=0;
particle_velocity_sum_smallest_half=0;
particle_velocity_sum_largest_half=0;
smallest_half = d<PARAMS.d0;
largest_half  = d>PARAMS.d0;
for it=1:Itlowres
    if t_vec_low_tres(it)>2
    particle_velocity_sum_all           = particle_velocity_sum_all           + sum(v_vec_low_tres(1,:,it));
    particle_velocity_sum_smallest_half = particle_velocity_sum_smallest_half + sum(v_vec_low_tres(1,smallest_half,it));
    particle_velocity_sum_largest_half  = particle_velocity_sum_largest_half  + sum(v_vec_low_tres(1,largest_half ,it));
    end
end
% particle_velocity_sum_all
% particle_velocity_sum_smallest_half
% particle_velocity_sum_largest_half
disp(['smallest half contributes for ',num2str(particle_velocity_sum_smallest_half/particle_velocity_sum_all*100,'%2.1f'),'%']);
disp(['largest  half contributes for ',num2str(particle_velocity_sum_largest_half/particle_velocity_sum_all*100,'%2.1f'),'%']);

fig=figure('units','centimeters','color',[1 1 1],'position',[2 2 8.5 8.5],'name','sediment flux vs time');
ax = subplot(1,1,1);
set(ax,'fontsize',fs);
plot(tvec,Sediment_flux_time_series);
xlabel('$t$ [s]','interpreter','latex','fontsize',fs);
ylabel('$Q$ [m$^2$/s]','interpreter','latex','fontsize',fs);
if FLAG_export_figs
        export_fig('TEST/Figures/sediment_transport_rate_vs_time.pdf','-nocrop');
        export_fig('TEST/Figures/sediment_transport_rate_vs_time.png','-nocrop');
end

%%
max_d=max(d);
min_d=min(d);
tail_length = 2;
fig=figure('color',[1 1 1],'units','centimeters');
ax = subplot(1,1,1);
figpos = get(fig,'position');
figpos([3,4]) = [8.5,10]; %17.5; % 17.5  / 8.5  [cm] (double/single column)
set(fig,'position',figpos);
it=1;
for p=1:P       
    alpha = (d(p)-min_d)/(max_d-min_d);
%                         colorfunction = alpha*[0 0 1]+(1-alpha)*[1 0 0]+[0 0.1 0];
    myplot(r_vec_low_tres(:,p,it),d(p),'color',my_colorfunction(alpha));
    if p==1
        hold on;
    end
    if plot_tails==1; plot_tail(squeeze(r_vec_low_tres(1,p,max(it-tail_length,1):it)),squeeze(r_vec_low_tres(2,p,max(it-tail_length,1):it)),squeeze(v_vec_low_tres(1,p,max(it-tail_length,1):it)),Lx,'color',[1 1 1]*0.6); end
    if r_vec_low_tres(1,p,it)-Lx+d(p)/2>0
        myplot(r_vec_low_tres(:,p,it)-[Lx;0;0],d(p),'color',my_colorfunction(alpha));
    end
    if r_vec_low_tres(1,p,it)+Lx-d(p)/2<Lx
        myplot(r_vec_low_tres(:,p,it)+[Lx;0;0],d(p),'color',my_colorfunction(alpha));
    end       
end
hold off;
view(0,90);
xlabel('x [m]');
ylabel('z [m]');
title(['t = ',num2str(t_vec_low_tres(it)),' [s]']);
cbar = colorbar;
% c.Label.String = titl;
% c.Limits = [minC 500];
cbar.Ticks = [0:0.0001:0.1];
limits = cbar.Limits;
caxis([min_d max_d]);
cbar.Limits = [min_d max_d];
title(cbar,'d [m]');
cInt = linspace(limits(1),limits(2),64);
alpha = linspace(0,1,100); cmap = zeros(100,3);
for i=1:numel(alpha)
    cmap(i,:) = my_colorfunction(alpha(i));
end
colormap(cmap);
% cmap = flipud(my_colorfunction(linspace(0,1,100)));
%         axis square;
xlim([0 Lx]);
plotymax = 0.2;
axis equal;
ylims=ylim;
ylim(ylims-ylims(1));
% ylim([0 plotymax]);

%%

max_d=max(d);
min_d=min(d);
if show_animation
tail_length = 2;
fig=figure('units','centimeters','color',[1 1 1],'position',[2 2 8.5 8.5],'name','blank');
for it=1:50:Itlowres
    subplot(1,2,1);
        for p=1:P       
            alpha = (d(p)-min_d)/(max_d-min_d);
%                         colorfunction = alpha*[0 0 1]+(1-alpha)*[1 0 0]+[0 0.1 0];
            myplot(r_vec_low_tres(:,p,it),d(p),'color',my_colorfunction(alpha));
            if p==1
                hold on;
            end
            if plot_tails==1; plot_tail(squeeze(r_vec_low_tres(1,p,max(it-tail_length,1):it)),squeeze(r_vec_low_tres(2,p,max(it-tail_length,1):it)),squeeze(v_vec_low_tres(1,p,max(it-tail_length,1):it)),Lx,'color',[1 1 1]*0.6); end
            if r_vec_low_tres(1,p,it)-Lx+d(p)/2>0
                myplot(r_vec_low_tres(:,p,it)-[Lx;0;0],d(p),'color',my_colorfunction(alpha));
            end
            if r_vec_low_tres(1,p,it)+Lx-d(p)/2<Lx
                myplot(r_vec_low_tres(:,p,it)+[Lx;0;0],d(p),'color',my_colorfunction(alpha));
            end       
        end
        hold off;
        view(0,90);
        xlabel('x');
        ylabel('z');
        title(['t = ',num2str(t_vec_low_tres(it)),' [s]']);
        xlim([0 Lx]);
        plotymax = 0.2;
        ylim([0 plotymax]);
    subplot(1,2,2);
        plot(u_low_tres(:,1),z,'color',[0 0 0.4]);
        hold on;
        plot(u_low_tres(:,it),z);
        plot(v_vec_low_tres(1,:,it),r_vec_low_tres(2,:,it),'k.');
        hold off;
        xlabel('u [m/s]');
        ylabel('z [m]');
        xlim([-1 maxu]);
        ylim([0 plotymax]);
        drawnow;
end
end

figure('units','centimeters','color',[1 1 1],'position',[2 2 8.5 8.5],'name','blank');
for i=1:2
    if i==1 
        it = 1;
    else
        it = numel(t_vec_low_tres);
    end
       
    subplot(1,2,i);
            for p=1:P       
                            alpha = (d(p)-min_d)/(max_d-min_d);
    %                         colorfunction = alpha*[0 0 1]+(1-alpha)*[1 0 0]+[0 0.1 0];
                            myplot(r_vec_low_tres(:,p,it),d(p),'color',my_colorfunction(alpha));
    %                         myplot(r_vec_low_tres(:,p,it),d(p),'color',color(:,mod(p-1,7)+1)); 
                            if p==1
                                hold on;
                            end
                            if r_vec_low_tres(1,p,it)-Lx+d(p)/2>0
                                myplot(r_vec_low_tres(:,p,it)-[Lx;0;0],d(p),'color',my_colorfunction(alpha));
                            end
                            if r_vec_low_tres(1,p,it)+Lx-d(p)/2<Lx
                                myplot(r_vec_low_tres(:,p,it)+[Lx;0;0],d(p),'color',my_colorfunction(alpha));
                            end       
            end
            hold off;
            view(0,90);
            xlabel('x');
            ylabel('z');
            title(['t = ',num2str(t_vec_low_tres(it)),' [s]']);
            xlim([0 Lx]);
            plotymax = 0.2;
            ylim([0 Lx]);
end

function s = myplot(x,d,varargin)
    theta = linspace(0,2*pi,10);
    X = d/2*cos(theta);
    Y = d/2*sin(theta);
    s=fill(x(1)+X,x(2)+Y,varargin{2},'LineStyle','none');
%     s = plot(x(1)+X,x(2)+Y,varargin{:});      
end

function [rgb_values] = my_colorfunction(alpha)
    if alpha<0
        rgb_values = [0 0 1];
    elseif alpha<0.5
        rgb_values = (1-2*alpha)*[0 0 1]+alpha*[0 1 0];
    elseif alpha>=0.5 && alpha<=1.0
        rgb_values = (1-1*alpha)*[0 1 0]+(alpha-0.5)*2*[1 0 0];
    else
        rgb_values = [1 0 0];
    end
end

function u=u_loglaw(z,PARAMS,Cplus)
    u_star=PARAMS.u_star;
    nu_f = PARAMS.nu_f;
    kappa = 0.41;
    u = u_star/kappa*log(z*u_star/nu_f)+u_star*Cplus;
end

function Qperbin = calc_Q_perbin(d,pfixed,Nbins_d,Bincenters,v_vec_low_tres,t_vec_low_tres,start_t,Lx)
    Qperbin = zeros(Nbins_d,1);
    Binwidth = Bincenters(2)-Bincenters(1);
    [~,P,~]=size(v_vec_low_tres);
    ind_larger_than_start_t = t_vec_low_tres>start_t;
    It_lowres=sum(ind_larger_than_start_t);
%     dt = t_vec_low_tres(2)-t_vec_low_tres(1);
    for p=pfixed+1:P
        ind = d(p)>=Bincenters-Binwidth/2 & d(p)<=Bincenters+Binwidth/2;
        Vp= 4/3*pi*(d(p)/2)^3;
        Qp = sum(v_vec_low_tres(1,p,ind_larger_than_start_t))/It_lowres/(Lx*d(p))*Vp;
        Qperbin(ind)=Qperbin(ind)+Qp;
    end
end

function p=plot_binned_data(Y,X,varargin)
x = zeros(numel(X)*2+2,1);
y = zeros(numel(X)*2+2,1);
binwidth=X(2)-X(1);
for i = 1:numel(X)
    x((i-1)*2+1+[1,2]) = X(i)+0.5*[-binwidth, binwidth];
    y((i-1)*2+1+[1,2]) = Y(i);
end
x(1)=X(1)-binwidth/2;
x(end)=X(end)+binwidth/2;
p=plot(y,x,varargin{:});
end

function Qperbin_dz = calc_Q_perbin_dz(d,pfixed,Nbins_d,Nbins_z,Bincenters_d,Bincenters_z,v_vec_low_tres,r_vec_low_tres,t_vec_low_tres,start_t,Lx)
    Qperbin_dz = zeros(Nbins_d,Nbins_z);
    % remove the time instances before t_start;
    v_vec_low_tres=v_vec_low_tres(:,:,t_vec_low_tres>start_t);
    r_vec_low_tres=r_vec_low_tres(:,:,t_vec_low_tres>start_t);
    t_vec_low_tres=t_vec_low_tres(t_vec_low_tres>start_t);
%     start_t
    It_lowres= numel(t_vec_low_tres);
    Binwidth_d = Bincenters_d(2)-Bincenters_d(1);
    Binwidth_z = Bincenters_z(2)-Bincenters_z(1);
    P=numel(d);
    progress = 0;
    textprogressbar('calculating outputs: ');
    for p=pfixed+1:P
        if mod(p,round(P/100))==0
            progress = progress +1;
            textprogressbar(progress);
%             disp(['progress calc Q per bin is: ',num2str(progress),'%']); 
        end
        ind_d = d(p)>=Bincenters_d-Binwidth_d/2 & d(p)<=Bincenters_d+Binwidth_d/2;
        Vp= 4/3*pi*(d(p)/2)^3;
        for nz = 1:Nbins_z
            ind_t = r_vec_low_tres(2,p,:)>=Bincenters_z(nz)-Binwidth_z/2 & r_vec_low_tres(2,p,:)<=Bincenters_z(nz)+Binwidth_z/2;
            Qp = sum(v_vec_low_tres(1,p,ind_t))/It_lowres/(Lx*d(p))*Vp;
            Qperbin_dz(ind_d,nz) = Qperbin_dz(ind_d,nz) + Qp; 
        end
%         Qp = sum(v_vec_low_tres(1,p,ind_larger_than_start_t))/It_lowres/(Lx*d(p))*Vp;
%         Qperbin(ind)=Qperbin(ind)+Qp;
    end
    textprogressbar('done');
end

function p = my_tile_pcolor(X,Y,Z,varargin)
    dx = X(2)-X(1);
    dy = Y(2)-Y(1);
    x = zeros(numel(X)+1,1);
    y = zeros(numel(Y)+1,1);
    z = zeros(numel(X)+1,numel(Y)+1);
    z(1:end-1,1:end-1)=Z;
    x(1:end-1)=X-dx/2; x(end)=X(end)+dx/2;
    y(1:end-1)=Y-dy/2; y(end)=Y(end)+dy/2;
    [xx,yy]=ndgrid(x,y);
    p= pcolor(xx,yy,z,varargin{:});
end

function textprogressbar(c)
    % This function creates a text progress bar. It should be called with a 
    % STRING argument to initialize and terminate. Otherwise the number correspoding 
    % to progress in % should be supplied.
    % INPUTS:   C   Either: Text string to initialize or terminate 
    %                       Percentage number to show progress 
    % OUTPUTS:  N/A
    % Example:  Please refer to demo_textprogressbar.m
    % Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
    % Version: 1.0
    % Changes tracker:  29.06.2010  - First version
    % Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/
    %% Initialization
    persistent strCR;           %   Carriage return pesistent variable
    % Vizualization parameters
    strPercentageLength = 10;   %   Length of percentage string (must be >5)
    strDotsMaximum      = 10;   %   The total number of dots in a progress bar
    %% Main 
    if isempty(strCR) && ~ischar(c)
        % Progress bar must be initialized with a string
        error('The text progress must be initialized with a string');
    elseif isempty(strCR) && ischar(c)
        % Progress bar - initialization
        fprintf('%s',c);
        strCR = -1;
    elseif ~isempty(strCR) && ischar(c)
        % Progress bar  - termination
        strCR = [];  
        fprintf([c '\n']);
    elseif isnumeric(c)
        % Progress bar - normal progress
        c = floor(c);
        percentageOut = [num2str(c) '%%'];
        percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
        nDots = floor(c/100*strDotsMaximum);
        dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
        strOut = [percentageOut dotOut];

        % Print it on the screen
        if strCR == -1
            % Don't do carriage return during first run
            fprintf(strOut);
        else
            % Do it during all the other runs
            fprintf([strCR strOut]);
        end

        % Update carriage return
        strCR = repmat('\b',1,length(strOut)-1);

    else
        % Any other unexpected input
        error('Unsupported argument type');
    end
end

function plot_tail(x,y,u,Lx,varargin)
    x2end = x(2:end);
    xdiff = x2end-x(1:end-1);
    ind = find(xdiff<0);
    if L2norm(u)<0.05
    elseif numel(ind)==0
        plot(x,y,varargin{:});
    elseif max(x)-min(x)<Lx/4
        plot(x,y,varargin{:});
    elseif L2norm(u)<0.1
    else
        if numel([x(1:ind);x(ind+1)+Lx])==numel(y(1:ind+1))
            plot([x(1:ind);x(ind+1)+Lx],y(1:ind+1),varargin{:});               %'color',[1 0 0]);
            plot([x(ind)-Lx; x(ind+1:numel(x))],y(ind:numel(x)),varargin{:});  %'color',[0 0 1]);
        end
    end
end
