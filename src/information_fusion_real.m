function X_hat=information_fusion_real(Y,Settings)

% Get number of time instances
[~,Settings.N]=size(Y);

% Allocate memory
X_hat=zeros(9,Settings.N);




% Generate accelerometer location matrix
Settings.R=Settings.alpha*Settings.acc_geometry;

% Generate the H matrix in the signal model.
H=zeros(3*(Settings.Ns+Settings.Nw),6);
for k=1:Settings.Ns
   H(3*(k-1)+1:3*k,1:3)=-skew_sym(Settings.R(:,k));
   H(3*(k-1)+1:3*k,4:6)=eye(3);  
end

% Generate the measurement error covariance
Qs=Settings.sigma_s^2*eye(3*Settings.Ns);
Qw=(Settings.sigma_w*pi/180)^2*eye(3*Settings.Nw);
Q=blkdiag(Qs,Qw);
Qinv=inv(Q);


for n=1:Settings.N
    
    % Check if any gyroscopes are saturated and generate the associated
    % indicator matrix to remove these entries from the signal model
    ind=[all((abs(Y(3*Settings.Ns+1:3:end,n))<Settings.gamma_w*pi/180));...
        all((abs(Y(3*Settings.Ns+2:3:end,n))<Settings.gamma_w*pi/180)); ...
        all((abs(Y(3*Settings.Ns+3:3:end,n))<Settings.gamma_w*pi/180))];
    
    ind=logical([ones(3*Settings.Ns,1); ...
        repmat(ind,Settings.Nw,1)]);
       
    % Calculate the the weighted least square matrix
    WLS=(H(ind,:)'*Qinv(ind,ind)*H(ind,:))\(H(ind,:)'*Qinv(ind,ind));
    
    % Calculate the weight matrix P
    P=Qinv(ind,ind)-Qinv(ind,ind)*H(ind,:)*WLS;
    
    % Initial estimate
    w0=mean(reshape(Y(3*Settings.Ns+1:end,n),3,Settings.Nw),2);
    ind_sat=(abs(w0)>=Settings.gamma_w*pi/180);
    if any(ind_sat);
% %     x=tensor_method_real(Y(1:3*Settings.Ns,n),Settings);
% %     w=sign(w0).*x(1:3);
    w=X_hat(1:3,n-1);%+0.001*X_hat(4:6,n-1);
    w0(ind_sat)=w(ind_sat);
    end
%     
    
    % Estimate the angular velocity
    %y=[Y(1:3*Settings.Ns,n+1); Y(3*Settings.Ns+1:end,n)];
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
    J=Grad(w,Settings);
    %J=gradg(Settings.R,w)';
    
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