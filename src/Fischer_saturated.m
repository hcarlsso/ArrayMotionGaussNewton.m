function J = Fischer_saturated(w,r,Qinv,gyro_limit)

    N = size(r,2);
    
    notsaturated = w<gyro_limit & w>-gyro_limit;
    notsaturated = [logical(ones(3*N,1)); repmat(notsaturated,N,1)];

    dg = gradg(r,w);
    H = Hmatrix(r);
    
    dg = dg(:,notsaturated);
    H = H(notsaturated,:);
    Qinv2 = Qinv(notsaturated,notsaturated);
    
    J = [dg*Qinv2*dg' dg*Qinv2*H;
         H'*Qinv2*dg' H'*Qinv2*H];
return

%     N = size(r,2);
% 
%     mdgx = dgx(w,r);
%     mdgy = dgy(w,r);
%     mdgz = dgz(w,r);
%     dffx = [repmat([1 0 0]',N,1); zeros(3*N,1)];
%     dffy = [repmat([0 1 0]',N,1); zeros(3*N,1)];
%     dffz = [repmat([0 0 1]',N,1); zeros(3*N,1)];
%     dfdwpx = zeros(2*3*N,1);
%     dfdwpy = zeros(2*3*N,1);
%     dfdwpz = zeros(2*3*N,1);
%     for j=1:N
%         dfdwpx(j*3-2:j*3) = [0 r(3,j) -r(2,j)]';
%         dfdwpy(j*3-2:j*3) = [-r(3,j) 0 r(1,j)]';
%         dfdwpz(j*3-2:j*3) = [r(2,j) -r(1,j) 0]';
%     end
%     
%     
%     
%     J11 = [mdgx'*Qinv*mdgx mdgx'*Qinv*mdgy mdgx'*Qinv*mdgz
%            mdgy'*Qinv*mdgx mdgy'*Qinv*mdgy mdgy'*Qinv*mdgz;
%            mdgz'*Qinv*mdgx mdgz'*Qinv*mdgy mdgz'*Qinv*mdgz];
%     J12 = [mdgx'*Qinv*dffx mdgx'*Qinv*dffy mdgx'*Qinv*dffz;
%            mdgy'*Qinv*dffx mdgy'*Qinv*dffy mdgy'*Qinv*dffz;
%            mdgz'*Qinv*dffx mdgz'*Qinv*dffy mdgz'*Qinv*dffz];
%     J13 = [mdgx'*Qinv*dfdwpx mdgx'*Qinv*dfdwpy mdgx'*Qinv*dfdwpz;
%            mdgy'*Qinv*dfdwpx mdgy'*Qinv*dfdwpy mdgy'*Qinv*dfdwpz;
%            mdgz'*Qinv*dfdwpx mdgz'*Qinv*dfdwpy mdgz'*Qinv*dfdwpz];
%     J22 = [dffx'*Qinv*dffx dffx'*Qinv*dffy dffx'*Qinv*dffz;
%            dffy'*Qinv*dffx dffy'*Qinv*dffy dffy'*Qinv*dffz;
%            dffz'*Qinv*dffx dffz'*Qinv*dffy dffz'*Qinv*dffz];
%     J23 = [dffx'*Qinv*dfdwpx dffx'*Qinv*dfdwpy dffx'*Qinv*dfdwpz;
%            dffy'*Qinv*dfdwpx dffy'*Qinv*dfdwpy dffy'*Qinv*dfdwpz;
%            dffz'*Qinv*dfdwpx dffz'*Qinv*dfdwpy dffz'*Qinv*dfdwpz];
%     J33 = [dfdwpx'*Qinv*dfdwpx dfdwpx'*Qinv*dfdwpy dfdwpx'*Qinv*dfdwpz;
%            dfdwpy'*Qinv*dfdwpx dfdwpy'*Qinv*dfdwpy dfdwpy'*Qinv*dfdwpz;
%            dfdwpz'*Qinv*dfdwpx dfdwpz'*Qinv*dfdwpy dfdwpz'*Qinv*dfdwpz];
%        
%     J = [J11 J12 J13;
%          J12' J22 J23;
%          J13' J23' J33];
% return