function [beta, nc] = calc_assoc_prob_tree(Omega, F, PDt, lambda, V)
%CALC_ASSOC_PROB_TREE Calculate table with all marginal probabilities of 
% association events for the joint probabilistic data association filter
% by an association tree
%
% Usage:
% [beta, nc] = calc_assoc_prob_tree(Omega, F, type, PDt, lambda, V)
%
% Inputs:
% Omega    = Validation matrix with all possible association events
% F        = Matrix of joint probabilities of measurements to targets
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

%% Tree-based association events
events = cell(mk,1);

% Indices of tracks
for j = 1:mk
    cond = (Omega(j,:) == 1);
    ind = find(cond);
    events{j,1} = ind;
end

%% JPDA with tree
% Create an empty tree
tr = tree('root');

% Initialize variables
n = 0;
n0 = 0;
k0 = 0;
nodes = zeros(c,1);
layers = cell(1,c);
edges = cell(c,c);

layers{1,1} = 0;

for j = 1:mk % For each layer - track - measurement
    for k = k0:n % Sweep all nodes from previous layer of measurements
        for i = 1:size(events{j,1},2) % For all valid tracks for each measurement
            if j == 1 % Addition of nodes from tree root
                [tr, nodes(i+1)] = tr.addnode(1, events{j,1}(1,i));
                % Number of added nodes from root
                n0 = n0+1;
                layers{1,n0+1} = j;
                edges{k+1,n0+1} = events{j,1}(1,i)-1;
            else % Addition of nodes from subsequent layer nodes
                % If target = 0 (no detection) or target was not used yet
                % by other measurement
                if events{j,1}(1,i) == 1 || ~sum(events{j,1}(1,i) == tr.Node{k+1,1})
                    [tr, nodes(n+2)] = tr.addnode(nodes(k+1,1), [tr.Node{k+1,1} events{j,1}(1,i)]);
                    % Number of added nodes from subsequent roots
                    n = n+1;
                    layers{1,n+1} = j;
                    edges{k+1,n+1} = events{j,1}(1,i)-1;
                end
            end
        end
    end
    if j > 1 % if not search from root
        k0 = n0+1; % start node indices from current layer
        n0 = n; % save current value of n for the next assignment to k0 = n0 + 1 = n + 1;
    else % if search from root
        k0 = 1; % start of node indices is 1
        n = n0; % end of node indices is the number of added nodes by first iteration
    end
end

% Total number of nodes
n = n+1;

% Get only the sequences described by the leafs (last n/2 nodes)
leaves = tr.findleaves;
k = length(leaves);

%tr.Layers = layers;
%tr.Edges = edges;
%figure('units','normalized','outerposition',[0 0 1 0.75],'Color','white');
%treeplot(tr.Parent',layers,edges);
%set(gcf, 'PaperPositionMode', 'auto', 'Color', 'w');
%set(gca,'LooseInset',get(gca,'TightInset'));
%export_fig('association_tree.pdf', '-pdf', '-native', '-nocrop');

% Allocate memory
Omf = zeros(k*mk,size(Omega,2));
deltaf = zeros(1,k*Nt);
tauf = zeros(k*mk,1);
phif = zeros(k,1);

% Enumerate all tables
for i = 1:k
    for j = 1:mk
        Omf(i*mk-mk+j,tr.Node{leaves(1,i),1}(1,j)) = 1;
    end
    tauf(i*mk -mk +1:i*mk,1) = sum(Omf(i*mk -mk +1:i*mk,2:end),2); % Defined for t = 1..Nt (exclude target 0)
    deltaf(1,i*Nt -Nt +1:i*Nt) = sum(Omf(i*mk -mk +1:i*mk,2:end),1); % Defined for t = 1..Nt (exclude target 0)
    phif(i,1) = sum((ones(mk,1)-tauf(i*mk -mk +1:i*mk,1)),1);
end
%% JPDA with tree -- end

beta = zeros(size(Omega));

% Target t = 0 (no detection) shall not be included
% Already taken into account in the probability of false detection
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

% Normalize marginals for each valid track
for j = 1:mk
    if sum(beta(j,:),2) > 0
        beta(j,:) = beta(j,:)/sum(beta(j,:),2);
    end
end

% Number of possible combinations
nc = k;

end

