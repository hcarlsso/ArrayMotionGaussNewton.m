function X_hat=tensor_method(Y,Settings,w)


% Allocate memory
X_hat=zeros(9,Settings.N);

% Generate accelerometer location matrix
R=[Settings.alpha*Settings.acc_geometry; ones(1,Settings.Ns)];

% Chaeck that 
if rank(R)~=4
   X_hat=NaN(9,Settings.N);
   return
end

% Generate projection matrix
P=R'/(R*R');



for n=1:Settings.N
    
    % Least square estimate
    T=reshape(Y(1:3*Settings.Ns,n),3,Settings.Ns)*P;
    
    % Angular velocity
    Ws=0.5*(T(1:3,1:3)+T(1:3,1:3)');
    
    X_hat(1:3,n)=sign(w(n)*Settings.w_norm).*...
        sqrt(abs(diag(Ws-0.5*trace(Ws)*eye(3))));
    
    % Angular acceleration
    X_hat(4,n)=0.5*(T(3,2)-T(2,3));
    X_hat(5,n)=0.5*(T(1,3)-T(3,1));
    X_hat(6,n)=0.5*(T(2,1)-T(1,2));
    
    % Specific force
    X_hat(7:9,n)=T(:,4);
end







end
