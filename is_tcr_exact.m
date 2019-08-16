% 'is_tcr_exact.m' measures the exactness of the tight-and-cheap relaxation
% for the ACOPF problem
% INPUTS
%   casedata: MATPOWER case
%   model: either 0 for loss minimization or 1 for cost minimization
function [optval, opttcr, error, optgap, solgap, isexact] = is_tcr_exact(casedata,model)
[Vtcr, vtcr, ~, ~, opttcr] = solve_opf_tcr(casedata,model);
vtcr = vtcr(:,1);
% Exactness error (epsilon)
error = 100*max(1 -  abs(vtcr)./sqrt(diag(Vtcr)));
[optval, optsol] = solve_opf(casedata,model);
% Optimality gap (gamma)
optgap = 100*(1 - opttcr/max(opttcr,optval));
% Optimality distance (rho)
solgap = 100*norm(optsol - vtcr)/norm(optsol);
isexact = (error < 0.01);