function X_hat=information_fusion(Y,Settings,w)
% But w is what we want to estimate
% Allocate memory
X_hat=zeros(9,Settings.N);


% -------- Same as generate_data.m ----------------

% Generate accelerometer location matrix
% Without any error 
% Settings is changed here, will it have an effect later on?
Settings.R=Settings.alpha*Settings.acc_geometry;

% Generate the H matrix in the signal model.
H=zeros(3*(Settings.Ns+Settings.Nw),6);
for k=1:Settings.Ns
   H(3*(k-1)+1:3*k,1:3)=-skew_sym(Settings.R(:,k));
   H(3*(k-1)+1:3*k,4:6)=eye(3);  
end

% Generate the measurement error covariance
% But Q is diagonal so the inverse is just Q_{ii}^{-1}
Qs=Settings.sigma_s^2*eye(3*Settings.Ns);
Qw=(Settings.sigma_w*pi/180)^2*eye(3*Settings.Nw);
Q=blkdiag(Qs,Qw);
Qinv=inv(Q);

% -------- Same as generate_data.m ----------------


% Iterate over time samples.
for n=1:Settings.N
    
    % Check if any gyroscopes are saturated and generate the associated
    % indicator matrix to remove these entries from the signal
    % model
    
    % Should the gyroscopes be removed? so size of vector should be
    % Settings.Nw?
    ind=logical([ones(3*Settings.Ns,1); ...
        ~(Y(3*Settings.Ns+1:end,n)>=Settings.gamma_w*pi/180)]);
       
    % Calculate the the weighted least square matrix
    WLS=(H(ind,:)'*Qinv(ind,ind)*H(ind,:))\(H(ind,:)'*Qinv(ind,ind));
    
    % Calculate the weight matrix P
    P=Qinv(ind,ind)-Qinv(ind,ind)*H(ind,:)*WLS;
    
    % Initial estimate
    % Time step n, in radians
    w0=pi/180.*w(:,n).*Settings.w_norm;
    
    % Estimate the angular velocity
    [X_hat(1:3,n), res]=GaussNewton(Y(:,n),w0,P,ind,Settings);
    
    % Estimate the angular acceleration and specific force
    X_hat(4:9,n)=WLS*res;
    
end


end


%% ------------------- Subfunctions -------------------------------------%%

function [w, res]=GaussNewton(y,w0,P,ind,Settings)

itr_ctr=0;
w=w0;
dw=inf(3,1);

while(itr_ctr<Settings.MaxIterations && norm(dw)>Settings.Tol)
    
    % Calculate the gradient and remove the entries of the saturated
    % gyroscopes
    J=Grad(w,Settings);%J=Gradient(w,Settings);
    J=J(ind,:);
    
    % Update the angular velocity estimate
    res=y-h(w,Settings);
    res=res(ind);
    dw=(J'*P*J)\(J'*P*res);
    w=w+dw;
    
    % Update the iteration counter
    itr_ctr=itr_ctr+1;
end

if itr_ctr==Settings.MaxIterations
    disp(itr_ctr)
end
end




function z=h(w,Settings)
z=[reshape((skew_sym(w)^2)*Settings.R,3*Settings.Ns,1);...
    kron(ones(Settings.Nw,1),w)];
end
