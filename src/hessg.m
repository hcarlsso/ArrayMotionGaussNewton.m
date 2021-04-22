function [Bx,By,Bz] = hessg(r)

N = size(r,2);

Bx = zeros(3,6*N);
ex = [1; 0; 0];
for i=1:N
    Bx(:,i*3-2:i*3) = r(:,i)'*ex*eye(3) + r(:,i)*ex' - 2*ex'*r(:,i);
end

By = zeros(3,6*N);
ey = [0; 1; 0];
for i=1:N
    By(:,i*3-2:i*3) = r(:,i)'*ey*eye(3) + r(:,i)*ey' - 2*ey'*r(:,i);
end

Bz = zeros(3,6*N);
ez = [0; 0; 1];
for i=1:N
    Bz(:,i*3-2:i*3) = r(:,i)'*ez*eye(3) + r(:,i)*ez' - 2*ez'*r(:,i);
end

ddg = [Bx By Bz];

end