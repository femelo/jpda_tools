function [betac, nc] = calc_assoc_prob_lbp(Omega, F, PDt, lambda, V)
%CALC_ASSOC_PROB_LBP Calculate table with all marginal probabilities of 
% association events for the joint probabilistic data association filter
% by loopy belief propagation
%
% Usage:
% [beta, nc] = calc_assoc_prob_lbp(Omega, F, PDt, lambda, V)
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
% University of Liverpool, February, 2014
%

if isempty(Omega) || isempty(F)
    betac = [];
    nc = 0;
    return;
end

%% Generation of validation matrices for all possible combinations of valid events

mk = size(Omega,1);
Nt = size(Omega,2)-1;

P = zeros(size(F));

% P(wjt)
for t = 1:Nt+1
    for j = 1:mk
        P1 = PDt*F(j,t)/lambda;
        P2 = (1-PDt);
        P(j,t) = (t > 1)*P1 + (t == 1)*P2;
    end
end

%% Jason Willians code

% Initialisation
w = P(:,2:end)'/(1-PDt);

n = mk;
m = Nt;
om = ones(1,m);
on = ones(1,n);
muba = ones(m,n);

muab = zeros(size(w));
muba0 = 0.1*muba;
muab0 = 0.1*ones(size(w));

while max(max(abs(muba0-muba))) > 1e-9 && max(max(abs(muab0-muab))) > 1e-9
    muba0 = muba;
    muab0 = muab;
    prodfact = muba .* w;
    sumprod = 1 + sum(prodfact);
    muab = w ./ (sumprod(om,:) - prodfact);
    summuab = 1 + sum(muab,2);
    muba = 1 ./ (summuab(:,on) - muab);
end

bf = [on; (muba .* w)];
betac = (bf./(ones(m+1,1)*sum(bf,1)))';

nc = 0;

end

