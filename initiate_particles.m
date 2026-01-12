function [r,v,omega] = initiate_particles(P,Lx,Ly,d)
    
    max_inner_itter = 100;
    max_outer_itter = 50;    
    good_solution = 0;
    outer_itter=1;
    while good_solution==0 && outer_itter<=max_outer_itter
        r      = zeros(3,P);
        v      = zeros(3,P);
        omega  = zeros(3,P);
        for p=1:P
            check=1;
            itter = 0;
            while check==1 && itter<max_inner_itter
                r(:,p)  = [mod(rand(1)*Lx,Lx); mod(rand(1)*Ly,0.5*Ly)+d(p)/2+d(1)/2; 0];
                max_delta = 0;
                    for i=1:p-1
                        [delta,~] = calc_delta(d(i),d(p),r(:,i),r(:,p),Lx);
                        max_delta=max(delta,max_delta);
                    end
                if max_delta<=0
                    check=0;
                end
                itter = itter+1;
            end  
            if itter==max_inner_itter
                disp(['Warning: (outer_itter = ',num2str(outer_itter),'/',num2str(max_outer_itter),') particle p=',num2str(p),' could not be placed']);
                break
            end
            v(:,p)     = [0;0;0];
            omega(:,p) = [0;0;0];
        end
        outer_itter=outer_itter+1;
        if p==P && itter<max_inner_itter
            good_solution=1;
        end
    end
    if outer_itter>=max_outer_itter && itter>=max_inner_itter
        figure;
        theta = linspace(0,2*pi,20);
        for i=1:p
            x=r(1,i)+d(i)/2*cos(theta);
            y=r(2,i)+d(i)/2*sin(theta);
            plot(x,y); 
            if i==1; hold on; end
        end
        xlim([0 Lx]);
        ylim([0 Ly]);
        error('initiate particles could not places. Reduce number of particles or increase domain size.');
    end
end