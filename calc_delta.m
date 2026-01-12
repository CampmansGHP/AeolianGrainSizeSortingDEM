function [delta,N] = calc_delta(d1,d2,r1,r2,L)
    if abs(r1(1)-r2(1))<L/2
    elseif r1(1)-r2(1)<0
        r2(1) = r2(1)-L;
    elseif r1(1)-r2(1)>0
        r2(1) = r2(1)+L;
    else
        disp('error, I did not anticipate for this! (calc_delta)');
        if isnan(r1(1))
            disp('the possition of the particle has a NaN value');
            error('the positions of the particle has a NaN value');
        end
    end
    N = (r1-r2)/norm(r1-r2,2);
    dr =(r1-r2)'*N;
    delta = 0.5*(d1+d2)-dr;
end