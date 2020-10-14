%% clear memory, screen, and close all figures
clear, clc, close all;
addpath('export_fig');
addpath('matlab-tree');

%% Constant velocity model parameters
T = 1;
F2 = [1 T; 0 1];
% G2 = [T^2/2; T];
% Q2 = G2*G2';
Q2 = [T^3/3 T^2/2; T^2/2 T];

Fkm1 = blkdiag(F2,F2);
Hk = [1 0 0 0; 0 0 1 0];

%% Process equation for targets x[k] = sys_f(x[k-1], P[k-1], u[k-1], Q[k-1]);
nt = 05; % number of targets
nx = 04; % number of states
sys = @state_eq2;


%% Observation equation for targets z[k] = obs_f(x[k], P[k], v[k], R[k]);
nz = 2; % number of observations
obs = @output_eq2;

%% PDF of process noise and noise generator function
nu = 4;
Q0 = 1e-2*blkdiag(Q2,Q2);
Q = 1e-3*blkdiag(Q2,Q2);
% Q0 = 1e-3*blkdiag(Q2,Q2);
% Q = 1e-4*blkdiag(Q2,Q2);
mux = zeros(nx,1);
p_sys_noise   = @(u, Qu) mvnpdf(u, mux, Qu);
gen_sys_noise = @(Qu) mvnrnd(mux, Qu)';

%% PDF of observation noise and noise generator function
nv = 2;
% R = 1e-3*eye(nv);
R = 0.5*eye(nv);
% Rf = 0.05*eye(nv);
Rf = R;
muz = zeros(nz,1);
p_obs_noise   = @(v, Rv) mvnpdf(v, muz, Rv);
gen_obs_noise = @(Rv) mvnrnd(muz, Rv)';

%% Initial state covariance matrix
% P0 = zeros(nx);
P0 = Q0;
S0 = Hk*(Fkm1*P0*Fkm1' + Q0)*Hk' + R;

%% Number of time steps
T = 50;

%% Separate memory space
sys_f = cell(1,nt);
obs_f = cell(1,nt);

x = cell(T,nt);
xt = cell(T,nt);
cov_x = cell(T,nt);
cov_sys = cell(1,nt);

z = cell(T,nt);
zt = cell(T,nt);
cov_z = cell(T,nt);
cov_obs = cell(1,nt);

u = zeros(nu,T);
v = zeros(nv,T);

ap = -5;
bp = +5;
av = -2;
bv = +2;

%% Attribute parameters and simulate system for all targets
for t = 1:nt
    % Assign the state and output functions to the target
    sys_f{t} = sys;
    obs_f{t} = obs;
    
    % Assign the process and measurement noise covariances to the target
    cov_sys{t} = Q;
    cov_obs{t} = R;
    
    % Initialize
    u(:,1) = gen_sys_noise(Q0);     % initial process noise
    v(:,1) = gen_obs_noise(R);      % initial observation noise
    
    x0 = zeros(nx,1);
    x0([1 3],1) = ap + (bp-ap).*rand(2,1); % initial position for target t
    x0([2 4],1) = av + (bv-av).*rand(2,1); % initial velocity for target t
    x{1,t} = x0 + gen_sys_noise(Q);
    z{1,t} = obs(x0, P0, v(:,1), R);
    
    % True state and observation
    xt{1,t} = x{1,t};
    zt{1,t} = obs(x{1,t}, zeros(nx,nx), zeros(nz,1), zeros(nz,nz));
    
    for k = 2:T
        u(:,k) = gen_sys_noise(Q);     % process noise
        if k < T
            v(:,k) = gen_obs_noise(R);
        else
            v(:,k) = gen_obs_noise(Rf);
        end
        x{k,t} = sys_f{t}(x{k-1,t}, zeros(nx,nx), u(:,k), zeros(nx,nx));   % simulate state
        % x{k,t} = sys_f{t}(x{k-1,t}, zeros(nx,nx), zeros(nx,1), zeros(nx,nx));   % simulate state
        z{k,t} = obs_f{t}(x{k,t}, zeros(nx,nx), v(:,k), zeros(nz,nz));     % simulate observation
        
        % True state and observation (without noise)
        % xt{k,t} = sys_f{t}(xt{k-1,t}, zeros(nx,nx), zeros(nx,1), zeros(nx,nx));
        xt{k,t} = x{k,t};
        zt{k,t} = obs_f{t}(x{k,t}, zeros(nx,nx), zeros(nz,1), zeros(nz,nz));
    end
    
end

%% Allocate memory
xh = cell(T,nt);
zh = cell(T,nt);

for t = 1:nt
    xh{1,t} = x{1,t};
    % zh{1,t} = obs_f{t}(x{1,t}, 0, 0, 0);
    zh{1, t} = z{1,t};
    cov_x{1,t} = P0;
    
    % xh{1,t} = [z{2,t}(1,1); (z{2,t}(1,1)-z{1,t}(1,1))/T; z{2,t}(2,1); (z{2,t}(2,1)-z{1,t}(2,1))/T];
    % cov_x{1,t} = blkdiag([R(1,1), R(1,1)/T; R(1,1)/T, 2*R(1,1)/T^2], [R(2,2), R(2,2)/T; R(2,2)/T, 2*R(2,2)/T^2]);
end

% P = blkdiag([R(1,1), R(1,1)/T; R(1,1)/T, 2*R(1,1)/T^2], [R(2,2), R(2,2)/T; R(2,2)/T, 2*R(2,2)/T^2]);
% S = Hk*(Fkm1*P*Fkm1' + Q0)*Hk' + R;
P = P0;
S = S0;

%% Parameters
% Volume of validation region
gamma_ = chi2inv(0.99,nz);
lambda = 0.01;
% cnz = pi^(nz/2)/gamma(nz/2 + 1);
% Estimation of the total surveillance region
% (union of validation region for all targets)
% Vt = nt*cnz*sqrt(det(gamma_*(S)));
% lambda = 1/Vt;

params.k            = 1;                % initial iteration number
params.m            = nt;               % number of tracks
params.Nt           = nt;               % number of targets
params.cov_sys      = cov_sys;          % process noise covariance matrix
params.cov_obs      = cov_obs;          % measurement noise covariance matrix
params.PDt          = 0.95;             % detection probability of target
params.lambda       = lambda;           % spatial density of false measurements / clutter density
params.gamma        = chi2inv(0.99,nz); % gate threshold - probability (PG) - for confidence of 99% and nz degrees of freedom

% Nr Monte Carlo runs
Nr = 1;
NEES = zeros(T,1);
ERMS = zeros(T,1);
tic
for i = 1:Nr
    
    fprintf('Run = %d/%d\n',i,Nr);
    
    % Estimate state
    for k = 2:T
        % fprintf('Iteration = %d/%d\n',k,T);
        
        % State estimation and filtered observation
        params.k = k;
        % [xh(k,:), cov_x(k,:), zh(k,:)] = jpda_filter(sys_f, obs_f, xh(k-1,:), cov_x(k-1,:), z(k-1,:), params, 'non-parametric');
        [xh(k,:), cov_x(k,:), zh(k,:)] = jpda_filter(sys_f, obs_f, xh(k-1,:), cov_x(k-1,:), z(k-1,:), params, 'parametric');
        % [xh(k,:), cov_x(k,:), zh(k,:)] = jpda_filter(sys_f, obs_f, xh(k-1,:), cov_x(k-1,:), z(k-1,:), params, 'tree');
        % [xh(k,:), cov_x(k,:), zh(k,:)] = jpda_filter(sys_f, obs_f, xh(k-1,:), cov_x(k-1,:), z(k-1,:), params, 'lbp');
        
        % Computation of NEES
        NEESkt = 0;
        ERMSkt = 0;
        for t = 1:nt
            % NEESkt = NEESkt ...
            %     + 0.5*(xt{k,t} - xh{k,t})'*(cov_x{k,t}\(xt{k,t} - xh{k,t})) -0.5*nx*nt;
            NEESkt = NEESkt ...
                + (xt{k,t} - xh{k,t})'*(cov_x{k,t}\(xt{k,t} - xh{k,t}))/nx;
            ERMSkt = ERMSkt ...
                + sum((zt{k,t} - zh{k,t}).^2);
        end
        NEES(k,1) = NEES(k,1) + NEESkt/nt;
        ERMS(k,1) = ERMS(k,1) + ERMSkt/nt;
    end
end
timerVal = tic;
toc

NEES = NEES/Nr;
ERMS = sqrt(ERMS/Nr);

xv = zeros(T,nx,nt);
zv = zeros(T,nz,nt);
xhv = zeros(T,nx,nt);
zhv = zeros(T,nz,nt);

for t = 1:nt
    for k = 1:T
    xv(k,1:nx,t) =  xt{k,t};
    zv(k,1:nz,t) =  zt{k,t};
    xhv(k,1:nx,t) =  xh{k,t};
    zhv(k,1:nz,t) =  zh{k,t};
    end
end

% Plot
figure(1)
hnd = zeros(nt,1);
color = rand(3,1);
plot(zv(:,1,1), zv(:,2,1), 'Color', color); hold on;
hnd(1) = plot(zhv(:,1,1), zhv(:,2,1), 'o', 'MarkerFaceColor', color);
lbl = cell(nt,1);
lbl{1,1} = 'Target 1';
for t = 2:nt
    color = rand(3,1);
    plot(zv(:,1,t), zv(:,2,t), 'Color', color);
    hnd(t) = plot(zhv(:,1,t), zhv(:,2,t), 'o', 'MarkerFaceColor', color);
    lbl{t,1} = sprintf('Target %d', t);
end
legend(hnd,lbl);
set(gca,'FontSize',12);
title('State-space','FontSize',14);
xlabel('Coordinate X (m)','FontSize',14);
ylabel('Coordinate Y (m)','FontSize',14);

figure(2)
plot((1:T)',NEES);
ylabel('NEES','FontSize',14);
xlabel('Epoch','FontSize',14);
set(gca,'FontSize',12);
title(sprintf('Normalised Estimation Error Squared (NEES) - %d targets',nt),'FontSize',14);
grid on;

figure(3)
plot((1:T)',ERMS);
ylabel('RMSE (m)','FontSize',14);
xlabel('Epoch','FontSize',14);
set(gca,'FontSize',12);
title(sprintf('Root-Mean-Square Error (RMSE) - %d targets',nt),'FontSize',14);
grid on;

return;
