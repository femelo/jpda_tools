%% clear memory, screen, and close all figures
clear, clc, close all;
addpath('export_fig');
addpath('tree');

%% Process equation for targets x[k] = sys_f(x[k-1], P[k-1], u[k-1], Q[k-1]);
nt = 10; % number of targets
nx = 2; % number of states
sys = @state_eq;


%% Observation equation for targets z[k] = obs_f(x[k], P[k], v[k], R[k]);
nz = 2; % number of observations
obs = @output_eq;

%% PDF of process noise and noise generator function
nu = 2;
Q0 = 1e-2*eye(nu);
Q = 1e-3*eye(nu);
mu = [0; 0];
p_sys_noise   = @(u, Qu) mvnpdf(u, mu, Qu);
gen_sys_noise = @(Qu) mvnrnd(mu, Qu)';

%% PDF of observation noise and noise generator function
nv = 2;
R = 1e-3*eye(nu);
Rf = 0.05*eye(nu);
p_obs_noise   = @(v, Rv) mvnpdf(v, mu, Rv);
gen_obs_noise = @(Rv) mvnrnd(mu, Rv)';

%% Initial state covariance matrix
P0 = zeros(2);

%% Number of time steps
T = 40;

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

a = -1;
b = +1;

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
    
    x0 = a + (b-a).*rand(2,1);      % initial state (position) for target t
    x{1,t} = x0;
    z{1,t} = obs(x0, P0, v(:,1), R);
    
    % True state and observation
    xt{1,t} = x0;
    zt{1,t} = obs(x0, zeros(nx,nx), zeros(nz,1), zeros(nz,nz));
    
    for k = 2:T
        u(:,k) = gen_sys_noise(Q);     % process noise
        if k < T
            v(:,k) = gen_obs_noise(R);
        else
            v(:,k) = gen_obs_noise(Rf);
        end
        x{k,t} = sys_f{t}(x{k-1,t}, zeros(nx,nx), u(:,k), zeros(nx,nx));   % simulate state
        z{k,t} = obs_f{t}(x{k,t}, zeros(nx,nx), v(:,k), zeros(nz,nz));     % simulate observation
        
        % True state and observation (without noise)
        xt{k,t} = sys_f{t}(xt{k-1,t}, zeros(nx,nx), zeros(nx,1), zeros(nx,nx));
        zt{k,t} = obs_f{t}(xt{k,t}, zeros(nx,nx), zeros(nz,1), zeros(nz,nz));
    end
    
end

%% Separate memory
xh = cell(T,nt);
zh = cell(T,nt);

for t = 1:nt
    xh{1,t} = x{1,t};
    % zh{1,t} = obs_f{t}(x{1,t}, 0, 0, 0);
    zh{1, t} = z{1,t};
    cov_x{1,t} = P0;
end

params.k            = 1;                % initial iteration number
params.m            = nt;               % number of tracks
params.Nt           = nt;               % number of targets
params.cov_sys      = cov_sys;          % process noise covariance matrix
params.cov_obs      = cov_obs;          % measurement noise covariance matrix
params.PDt          = 0.8;              % detection probability of target
params.lambda       = 0.3317;           % spatial density of false measurements / clutter density
params.gamma        = chi2inv(0.99,nz); % gate threshold - probability (PG) - for confidence of 99% and nz degrees of freedom

%% Estimate state
for k = 2:T
   fprintf('Iteration = %d/%d\n',k,T);
   
   % State estimation and filtered observation
   params.k = k;
   % [xh(k,:), cov_x(k,:), zh(k,:)] = jpda_filter(sys_f, obs_f, xh(k-1,:), cov_x(k-1,:), z(k-1,:), params, 'parametric');
   [xh(k,:), cov_x(k,:), zh(k,:)] = jpda_filter(sys_f, obs_f, xh(k-1,:), cov_x(k-1,:), z(k-1,:), params, 'tree');  
 
end

xv = zeros(T,nx,nt);
zv = zeros(T,nz,nt);
xhv = zeros(T,nx,nt);
zhv = zeros(T,nz,nt);

for t = 1:nt
    for k = 1:T
%     xv(k,1:nx,t) =  xt{k,t};
%     zv(k,1:nz,t) =  zt{k,t};
    xv(k,1:nx,t) =  x{k,t};
    zv(k,1:nz,t) =  z{k,t};
    xhv(k,1:nx,t) =  xh{k,t};
    zhv(k,1:nz,t) =  zh{k,t};
    end
end

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

return;
