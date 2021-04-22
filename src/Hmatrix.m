
function H = Hmatrix(r)

N = size(r,2);

H = zeros(2*3*N,6);

for i=1:N
    H(i*3-2:i*3,:) = [eye(3) -skew(r(:,i))];
end

return