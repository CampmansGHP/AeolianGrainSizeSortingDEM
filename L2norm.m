function [L2]=L2norm(u)
    I=length(u(:));
    L2=sqrt(sum(u(:).*conj(u(:)))/I);
end