
function dg = gradg(r,w)

N = size(r,2);

dg = zeros(3,6*N);
for i=1:N
    dg(:,i*3-2:i*3) = w'*r(:,i)*eye(3)+r(:,i)*w'-2*w*r(:,i)';
    dg(:,i*3-2+3*N:i*3+3*N) = eye(3);
end

end
