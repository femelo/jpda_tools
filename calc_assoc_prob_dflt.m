function [beta, nc] = calc_assoc_prob_dflt(Omega, F, type, PDt, lambda, V)
%CALC_ASSOC_PROB_DFLT Calculate table with all marginal probabilities of 
% association events for the joint probabilistic data association filter
% by default method
%
% Usage:
% [beta, nc] = calc_assoc_prob_dflt(Omega, F, type, PDt, lambda, V)
%
% Inputs:
% Omega    = Validation matrix with all possible association events
% F        = Matrix of joint probabilities of measurements to targets
% type     = method of calculation of the marginal probabilities for assignments:
%   'parametric'     = default JPDAF with parametric clutter model
%   'non_parametric' = default JPDAF with parametric clutter model
% PDt      = detection probability of target
% lambda   = spatial density of false measurements / clutter density
% V        = volume of validation region (row vector - for each target)
%
% Outputs:
% beta     = Matrix of marginal probabilities of association events
% nc       = Number of combinations for valid events
%
% Coded by:
% Flavio Eler de Melo (flavio.eler@gmail.com)
% University of Liverpool, August, 2013
%

%% Generation of validation matrices for all possible combinations of valid events

mk = size(Omega,1);
Nt = size(Omega,2)-1;

% Combinations of unique (feasible) associations for all tracks (all rows)
sc = sum(Omega,2);
% All possible combinations without the constraint of at most one
% measurement originaged from a target
c = prod(sc);

% Association events
% theta = cell(c,1);

% Generation of feasible events
ind = cell(mk,1);
Om = zeros(c*mk,size(Omega,2));
Oms = zeros(c*mk,size(Omega,2));

delta = zeros(1,c*Nt);
tau = zeros(c*mk,1);
phi = zeros(c,1);

% Indices of tracks
for j = 1:mk
    ind{j} = find(Omega(j,:) == 1);
end

% Matrix with all association events that satisfy:
% one source for each measurement (one target per row)
for j = 1:mk % 1 to 4
    nt = sc(j); % 4 | 3 | 3 | 2
    nl = prod(sc(1:j)); % 4 | 4*3 (=12) | 4*3*3 (=36) | 4*3*3*2 (=72)
    pf = c*mk/nl;     % 72c*4m/4t (=72) | 72c*4m/(4*3)t (=24) 
    % 72c*4m/(4*3*3)t (=8) | 72c*4m/(4*3*3*2)t (=4)
    base = repmat(eye(nt),nl/nt,1); % eye(4) | eye(3) x 4 | eye(3) x 4*3 | eye(2) x 4*3*3
    for k = 1:nl % 1 to 4 | 1 to 4*3 | 1 to 4*3*3 | 1 to 4*3*3*2
        for L = 1:pf/mk
            Om(pf*k -pf +j +mk*L -mk,ind{j}) = base(k,1:nt);
        end
    end
end

% Eliminate combinations with more than one
% measurement originaged from a target
k = 0;
for i = 1:c
    Omi = Om(i*mk -mk +1:i*mk,:);
    deltai = sum(Omi(:,2:end),1); % Defined for t = 1..Nt (exclude dummy target 0)
    if ~(deltai > 1)
        k = k+1;
        Oms(k*mk -mk +1:k*mk,:) = Omi;
        tau(k*mk -mk +1:k*mk,1) = sum(Omi(:,2:end),2); % Defined for t = 1..Nt (exclude dummy target 0)
        delta(1,k*Nt -Nt +1:k*Nt) = deltai;
        phi(k,1) = sum((ones(mk,1)-tau(k*mk -mk +1:k*mk,1)),1);
        % theta{k,1} = Omi;
    end
end

% Allocate memory
Omf = zeros(k*mk,size(Omega,2));
deltaf = zeros(1,k*Nt);
tauf = zeros(k*mk,1);
phif = zeros(k,1);
% thetaf = cell(k,1);

% Transfer only the feasible events and parameters
Omf(1:k*mk,:) = Oms(1:k*mk,:);
deltaf(1,1:k*Nt) = delta(1,1:k*Nt);
tauf(1:k*mk,1) = tau(1:k*mk,1);
phif(1:k,1) = phi(1:k,1);
% thetaf(1:k,1) = theta(1:k,1);

beta = zeros(size(Omega));

% Target t = 0 (no detection) shall not be included
% Already taken into account in the probability of false detection
switch lower(type)
    case 'parametric'
        for t = 1:Nt+1
            for j = 1:mk
                betatj = 0;
                % ctj = 0;
                for i = 1:k
                    P1 = 1; P2 = 1;
                    % Only perform the calculation if it is part of the event to
                    % compute marginal probability beta(j,t)
                    if Omf(i*mk -mk +j,t) == 1
                        % For targets where
                        ind = find(sum(Omf(i*mk -mk +1:i*mk,:),1) == 1);
                        for ti = ind
                            a1 = ((lambda^-1)*F(:,ti)).^(tauf(i*mk -mk +1:i*mk,1).*Omf(i*mk -mk +1:i*mk,ti));
                            P1 = P1*prod(a1);
                        end
                        
                        a2 = (PDt.^deltaf(1,i*Nt -Nt +1:i*Nt)).*...
                            ((1 -PDt).^(1 -deltaf(1,i*Nt -Nt +1:i*Nt)));
                        P2 = prod(a2);
                        
                        P = P1*P2;
                        % Include if the assignment is such that wjt = 1
                        % if Omf(i*mk -mk +j,t) == 1
                        % if Omega(j,t) == 1
                        betatj = betatj + P;
                        % end
                        % Normalization (sum over all valid assignments)
                        % ctj = ctj + P;
                    end
                end
                beta(j,t) = betatj;
            end
        end
    case 'non-parametric'
        % Total volume of the surveillance region
        Vt = sum(V);
        for t = 1:Nt+1
            for j = 1:mk
                betatj = 0;
                % ctj = 0;
                for i = 1:k
                    P1 = 1; P2 = 1;
                    % Just perform the calculation if it is part of the event to
                    % compute marginal
                    if Omf(i*mk -mk +j,t) == 1
                        % For targets where
                        ind = find(sum(Omf(i*mk -mk +1:i*mk,:),1) == 1);
                        for ti = ind
                            a1 = (Vt*F(:,ti)).^(tauf(i*mk -mk +1:i*mk,1).*Omf(i*mk -mk +1:i*mk,ti));
                            P1 = P1*prod(a1);
                        end
                        
                        a2 = (PDt.^deltaf(1,i*Nt -Nt +1:i*Nt)).*...
                            ((1 -PDt).^(1 -deltaf(1,i*Nt -Nt +1:i*Nt)));
                        P2 = prod(a2);
                        
                        P = factorial(phif(i,1))*P1*P2;
                        
                        % Include if the assignment is such that wjt = 1
                        % if Omf(i*mk -mk +j,t) == 1
                        % if Omega(j,t) == 1
                        betatj = betatj + P;
                        % end
                        % Normalization (sum over all valid assignments)
                        % ctj = ctj + P;
                    end
                end
                beta(j,t) = betatj;
            end
        end
    otherwise
        error('Error: Unknown type.');
end

% Normalize marginals for each valid track
for j = 1:mk
    if sum(beta(j,:),2) > 0
        beta(j,:) = beta(j,:)/sum(beta(j,:),2);
    end
end

% Number of possible combinations
nc = k;

end

