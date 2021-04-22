%% Settings

% Name of the simulation case
name={'inplane','outofplane','all','allcube'};

% Direction of angular velocity vector
% Why are there four rows?
% 4x3 -> 3x4
w_direction=[1 0 0; 0 0 1; 1 1 1; 1 1 1]';

% Number of Monte Carlo loops
Settings.M=1e3;

% Angular velocity values
Settings.wmin=100; %[deg]
Settings.wmax=10000; %[deg]
Settings.N=100; %[Number of sample points]

% Standard deviation of sensor errors
Settings.sigma_s=0.01; % [m/s^2]
Settings.sigma_w=1;  % [deg/s]

% Saturation level gyroscopes
Settings.gamma_w=2000; % [deg]

% Geomentry
Settings.geometry=[];%'real';

% Standard deviation of the error in sensor location 
% Should maybe be variance?
Settings.s2r=0;%(1e-4)^2;  %[m]


% Gauss Newton settings
Settings.MaxIterations=30;
Settings.Tol=1e-6;


%% Loop over the different rotation directions and array geometries 
for k=2:2:4

    % Set the geometry of the array
    if strcmp(Settings.geometry,'real')
        load exp_data
        Settings.acc_geometry=PosMIMU4444BT(32);
        %ind=(Settings.acc_geometry(3,:)<0);
        %Settings.acc_geometry=Settings.acc_geometry(:,ind); 
        Settings.alpha=1;
        Settings.Ns=32;
        Settings.Nw=32;
        Settings.w_norm=w_unit; % What is this?
    else
        if k==4
            % 6 sensors?
            Settings.acc_geometry=[...
                -1 0 0; ...
                1 0 0; ...
                0 -1 0; ...
                0 1 0; ...
                0 0 1; ...
                0 0 -1]'; % 6x3 -> 3x6
            Settings.Ns=size(Settings.acc_geometry,2);
            Settings.Nw=Settings.Ns;
        else
            % 4 acc sensors?
            Settings.acc_geometry=[...
                -1 0 0; ...
                1 0 0; ...
                0 -1 0; ...
                0 1 0]';
            Settings.Ns=size(Settings.acc_geometry,2);
            Settings.Nw=Settings.Ns;
        end
        
        % Set the scale
        % scale factor for the positions?
        Settings.alpha=0.01;
        
        % Set the normalized angular velocity vector
        % Pick the k direction in w_direction
        Settings.w_norm=w_direction(:,k)./norm(w_direction(:,k));

    end



%% Monte Carlo loop

% Allocate memory
rmse_mle=zeros(9,Settings.N);
rmse_norm_mle=zeros(1,Settings.N);
rmse_tensor=zeros(9,Settings.N);

tic
for m=1:Settings.M
   
    % Display the iteration number
    if m/10==round(m/10)
        disp(m)
        toc
        tic
    end
    
    % Generate inertial sensor array data. The specific force and angular
    % acceleration is picked randomly as the performance of the estimator
    % should be indendent of these quantities.
    Settings.Angular_Acceleration=randn(3,1);
    Settings.Specific_Force=randn(3,1);
    [Y,w]=generate_data(Settings);
    
    % Run the proposed information fusion method
    x_hat_mle=information_fusion(Y,Settings,w);
    
    % Run the tensor based estimation method
    x_hat_tensor=tensor_method(Y,Settings,w);
    
    % Save the statistics
    rmse_mle(1:3,:)=rmse_mle(1:3,:)+(x_hat_mle(1:3,:)-pi/180.*Settings.w_norm*w).^2;
    rmse_mle(4:6,:)=rmse_mle(4:6,:)+(x_hat_mle(4:6,:)-Settings.Angular_Acceleration*ones(1,Settings.N)).^2;
    rmse_mle(7:9,:)=rmse_mle(7:9,:)+(x_hat_mle(7:9,:)-Settings.Specific_Force*ones(1,Settings.N)).^2;
    
    rmse_norm_mle=rmse_norm_mle+(sqrt(sum(x_hat_mle(1:3,:).^2))-pi/180.*w).^2;
    
    rmse_tensor(1:3,:)=rmse_tensor(1:3,:)+(x_hat_tensor(1:3,:)-pi/180.*Settings.w_norm*w).^2;
    rmse_tensor(4:6,:)=rmse_tensor(4:6,:)+(x_hat_tensor(4:6,:)-Settings.Angular_Acceleration*ones(1,Settings.N)).^2;
    rmse_tensor(7:9,:)=rmse_tensor(7:9,:)+(x_hat_tensor(7:9,:)-Settings.Specific_Force*ones(1,Settings.N)).^2;
end

% Calculate the rmse
rmse_mle=sqrt(rmse_mle./Settings.M);
rmse_tensor=sqrt(rmse_tensor./Settings.M);
rmse_norm_mle=sqrt(rmse_norm_mle./Settings.M);


%% Cramer-Rao bound

% CRB calculated using the general expression
Settings.gamma_w=inf;
CRB=CRB_general_case(Settings,w);

Settings.gamma_w=1;
CRB_sat=CRB_general_case(Settings,w);
Settings.gamma_w=2000;

%% Plot the results

figure(1)
clf
legend_h=zeros(1,8);
semilogy([Settings.gamma_w/Settings.w_norm(1) Settings.gamma_w/Settings.w_norm(1)],[10^(-1.3) 1e0],'k--','LineWidth',2)
hold on
semilogy([Settings.gamma_w/Settings.w_norm(2) Settings.gamma_w/Settings.w_norm(2)],[10^(-1.3) 1e0],'k--','LineWidth',2)
semilogy([Settings.gamma_w/Settings.w_norm(3) Settings.gamma_w/Settings.w_norm(3)],[10^(-1.3) 1e0],'k--','LineWidth',2)
semilogy([w(1) w(end)],Settings.sigma_w/sqrt(Settings.Nw)*ones(1,2),'k--','LineWidth',2)
legend_h(1:3)=semilogy(w,180/pi.*rmse_mle(1:3,:)');
axis([Settings.wmin Settings.wmax 10^(-1.3) 1e1])
legend_h(4)=semilogy(w,180/pi.*sqrt(CRB(1,:)),'k');
semilogy(w,180/pi.*sqrt(CRB(2,:)),'k')
semilogy(w,180/pi.*sqrt(CRB(3,:)),'k')
legend_h(5)=semilogy(w,180/pi.*sqrt(CRB_sat(1,:)),'k-.');
semilogy(w,180/pi.*sqrt(CRB_sat(2,:)),'k-.')
semilogy(w,180/pi.*sqrt(CRB_sat(3,:)),'k-.')
legend_h(6:8)=semilogy(w,180/pi.*rmse_tensor(1:3,:)','--');
grid on
box on
legend(legend_h,'x','y','z','CRB','CRB saturated gyroscope','x tensor','y tensor','z tensor')

text(7000,1.1*Settings.sigma_w/sqrt(Settings.Nw),'base level')
text(1.1*Settings.gamma_w,1.5*Settings.sigma_w/sqrt(Settings.Nw),'Gyroscope saturation level')
title('title')
xlabel('xaxis')
ylabel('yaxis')
title('Angular velocity')
xlabel('omega [deg]')
ylabel('rmse [deg/s]')


figure(2)
clf
semilogy(w,180/pi.*rmse_mle(4:6,:)')
axis([Settings.wmin Settings.wmax 1e1 1e2])
hold on
semilogy(w,180/pi.*rmse_tensor(4:6,:)','--');
semilogy(w,180/pi.*sqrt(CRB(4,:)),'k')
semilogy(w,180/pi.*sqrt(CRB(5,:)),'k')
semilogy(w,180/pi.*sqrt(CRB(6,:)),'k')
semilogy(w,180/pi.*sqrt(CRB_sat(4,:)),'k-.')
semilogy(w,180/pi.*sqrt(CRB_sat(5,:)),'k-.')
semilogy(w,180/pi.*sqrt(CRB_sat(6,:)),'k-.')
grid on
legend('x','y','z','x tensor','y tensor','z tensor','crb')

title('Angular acceleration')
xlabel('omega [deg/s]')
ylabel('rmse [deg/s^2]')




%% Save the result
%save(name{k})

end

