function [xhk, cov_xhk, zhk] = jpda_filter(sys_f, obs_f, xkm1, cov_xkm1, zk, params, type)
%% Parametric/Non-parametric joint probabilistic data association filter
%
% Usage:
% [xhk, cov_xhk, zhk] = jpda_filter(sys_f, obs_f, xkm1, cov_xkm1, zk, params, type)
%
%
% Inputs:
% sys_f    = cell array of function handles to state equations for each target (process)
% obs_f    = cell array of function handles to output equations for each target (measurement)
% xkm1     = cell array of state vectors at time k-1 (column vectors)
% cov_xkm1 = cell array of state covariance matrices at time k-1 
% zk       = cell array of observation vectors at time k (column vectors)
% params   = structure with the following fields
%   .k       = iteration number
%   .m       = number of measurements (tracks)
%   .Nt      = number of targets (known)
%   .cov_sys = process noise covariance matrix (Q for all targets)
%   .cov_obs = measurement noise covariance matrix (R for all targets)
%   .PDt     = detection probability of target
%   .lambda  = spatial density of false measurements / clutter density
%   .gamma   = gate threshold - probability (PG)
% type     = type of prior probability mass function of the number of false measurements (clutter model).
%            Set it either to 'parametric' or 'non_parametric'.
%
% Outputs:
% xhk     = cell array of estimated state vectors at time k (column vectors)
% cov_xhk = cell array of estimated state covariance matrices at time k
% zhk     = cell array of updated observations at time k
%
% Reference:
% [1] Bar-Shalom, Y., Willet, P.K., and Tian, X. Tracking and Data Fusion:
%     A Handbook of Algorithms. April 2011. YBS Publishing.
%     Chapters 3 and 6, p 174--202 and 385--405.
%
% Coded by:
% Flavio Eler de Melo (flavio.eler@gmail.com)
% University of Liverpool, August, 2013

%%
k = params.k;
if k == 1
   error('error: k must be an integer greater or equal than 2');
end

%% Initialize variables
m = size(zk,2);                          % number of measurements

Nt = params.Nt;                          % number of targets
cov_sys = params.cov_sys;                % process noise covariance matrix
cov_obs = params.cov_obs;                % measurement noise covariance matrix
PDt = params.PDt;                        % detection probability of target
lambda = params.lambda;                  % spatial density of false measurements / clutter density
gamma_ = params.gamma;                   % gate threshold probability (PG)

% Allocate memory
xhkp = cell(1,Nt);
cov_xhkp = cell(1,Nt);
zhkp = cell(1,Nt);
cov_zhkp = cell(1,Nt);
k_gain = cell(1,Nt);

xhk = cell(1,Nt);
cov_xhk = cell(1,Nt);
zhk = cell(1,Nt);

% Innovation vectors
nujt = cell(m,Nt);

% Dimension of state vector
nx = size(xkm1{1,1},1);
% Dimension of output vector
nz = size(zk{1,1},1);

% Volume of validation region
V = zeros(1,Nt+1);

% Initialize a maximum validation matrix
Omega = zeros(m,Nt+1);
Omega(:,1) = ones(m,1);

%% Prediction steps, calculation of innovations and validation
for t = 1:Nt
    % Prediction:
    % xh[k|k-1] = f( x[k|k-1] );
    % Ph[k|k-1] = F[k-1] . P[k-1|k-1] . F[k-1]' + Q[k-1]
    [xhkp{1,t}, cov_xhkp{1,t}] = sys_f{t}(xkm1{1,t}, cov_xkm1{1,t}, zeros(nx,1), cov_sys{1,t});
    % zh[k|k-1] = h( xh[k|k-1] )
    % S[k] = H[k] . Ph[k|k-1] . H[k]' + R[k]
    [zhkp{1,t}, cov_zhkp{1,t}, k_gain{1,t}] = obs_f{t}(xhkp{1,t}, cov_xhkp{1,t}, zeros(nz,1), cov_obs{1,t});
    
    % Validation
    for j = 1:m
        % V[k] = (z - zh[k|k-1])' . S[k]^-1 . (z - zh[k|k-1])
        Nuk = (zk{j} - zhkp{t})'*(cov_zhkp{t}\(zk{j} - zhkp{t}));
        % Validation
        if Nuk <= gamma_
            Omega(j,t+1) = 1;
        end
        nujt{j,t} = zk{j} - zhkp{t};
    end
    
    % Volume of validation region
    cnz = pi^(nz/2)/gamma(nz/2 + 1);
    V(1,t+1) = cnz*sqrt(det(gamma_*cov_zhkp{t}));
end

%% Number of validated measurements (sum over targets - columns)
indvm = find(sum(Omega,2) > 1)';
mk = length(indvm);
Omegaf = zeros(mk,Nt+1);
Omegaf(1:mk,:) = Omega(indvm,:);

% mk = m;
% Omegaf = zeros(mk,Nt+1);
% Omegaf(1:mk,:) = Omega(1:mk,:);

zkv = cell(mk,1);
zkv(1:mk) = zk(indvm);

nujtv = cell(mk,Nt);
nujtv(1:mk,:) = nujt(indvm,:);

%% Calculate joint association probabilities for targets/tracks
F = zeros(size(Omegaf));
F(:,1) = ones(size(F(:,1)));

for t = 2:Nt+1
	for j = 1:mk
        if Omegaf(j,t) == 1
            % F(j,t) = mvncdf(zkv{j}, zhkp{t-1}, cov_zhkp{t-1});
            % Try to avoid numerical asymmetry
            cov_z = (cov_zhkp{t-1} + cov_zhkp{t-1}.')/2;
            try
                F(j,t) = mvnpdf(zkv{j}, zhkp{t-1}, cov_z);
            catch err
                % Round covariance to 9 digits
                nd = 9; % number of digits / round covariance to 9 digits
                cov_z = round(cov_zhkp{t-1}*10^nd)./10^nd;
                F(j,t) = mvnpdf(zkv{j}, zhkp{t-1}, cov_z);
            end
        end
	end
end

[beta, k] = calc_assoc_prob(Omegaf, F, type, PDt, lambda, V);

%% For all targets, update state and state covariance
for t = 1:Nt
    nut = 0;
    nu2t = 0;
    
    if mk > 0    
        for j = 1:mk
            % Combined innovation
            % nu[k] = Sum_j{ betaj[k] . nuj[k] }
            nut = nut + beta(j,t+1)*nujtv{j,t};
            
            % Square weighted combined innovation
            % nu2[k] = Sum_j{ betaj[k] . nuj[k] . nuj[k]' }
            nu2t = nu2t + beta(j,t+1)*nujtv{j,t}*nujtv{j,t}';
        end
        
        % Update state
        % xh[k|k] = xh[k|k-1] + W[k] . nu[k];
        xhk{1,t} = xhkp{1,t} + k_gain{1,t}*nut;
        
        % Update state covariance matrix for correct measurement
        % Pc[k|k] = Ph[k|k-1] - W[k] . S[k] . W[k]'
        Pc = cov_xhkp{1,t} -k_gain{1,t}*cov_zhkp{1,t}*k_gain{1,t}';
        
        % Calculate spread of the innovations
        % Ps[k] = W[k] . ( nu2[k] - nu[k] . nu[k]' ) . W[k]'
        P = k_gain{1,t}*(nu2t - (nut*nut'))*k_gain{1,t}';
        
        % Update associated state covariance matrix
        % Ph[k|k] = beta0 . Ph[k|k-1] +(1 -beta0) . Pc[k|k] + Ps
        
        % Probability of none of measurements is originated from a target
        beta0 = 1 - sum(beta(1:mk,t+1));
        cov_xhk{1,t} = beta0*cov_xhkp{1,t} + (1 -beta0)*Pc + P;
        
    else
        % Update state - no valid innovation
        % xh[k|k] = xh[k|k-1] + W[k] . nu[k];
        xhk{1,t} = xhkp{1,t};
        
        % Update associated state covariance matrix - no valid innovation
        cov_xhk{1,t} = cov_xhkp{1,t};
    end
    
    % Update estimated measurement
    % zh[k|k] = h( xh[k|k] )
    try
        [zhk{1,t},~,~] = obs_f{t}(xhk{1,t}, cov_xhk{1,t}, zeros(nz,1), cov_obs{1,t});
    catch err
        error('Unknown error.');
    end
end

return
