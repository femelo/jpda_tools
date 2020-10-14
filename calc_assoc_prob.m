function [beta, nc] = calc_assoc_prob(Omega, F, type, PDt, lambda, V)
%CALC_ASSOC_PROB Calculate table with all marginal probabilities of 
% association events for the joint probabilistic data association filter
%
% Usage:
% [beta, nc] = calc_assoc_prob(Omega, F, type, PDt, lambda, V)
%
% Inputs:
% Omega    = Validation matrix with all possible association events
% F        = Matrix of joint probabilities of measurements to targets
% type     = method of calculation of the marginal probabilities for assignments:
%   'parametric'     = default JPDAF with parametric clutter model
%   'non-parametric' = default JPDAF with parametric clutter model
%   'tree'           = JPDAF with association tree
%   'lbp'            = JPDAF with loopy belief propagation
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

% Target t = 0 (no detection) shall not be included
% Already taken into account in the probability of false detection
switch lower(type)
    case {'parametric' 'non-parametric'}
        [beta, nc] = calc_assoc_prob_dflt(Omega, F, type, PDt, lambda, V);
    case 'tree'
        [beta, nc] = calc_assoc_prob_tree(Omega, F, PDt, lambda, V);
    case 'lbp'
        [beta, nc] = calc_assoc_prob_lbp(Omega, F, PDt, lambda, V);
    otherwise
        error('Error: Unknown type.');
end

end

