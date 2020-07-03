% 'is_tcr_qcr_exact.m' measures the exactness of TCR and QCR
% for the ACOPF problem
% INPUTS
%   casedata: MATPOWER case
%   model: either 0 for loss minimization or 1 for cost minimization
function [optgap, error, solgap, isexact] = is_tcr_qcr_exact(casedata,model)
[opttcr, soltcr, Vtcr, ~, ~] = solve_opf_tcr(casedata,model);
[optqcr, solqcr, Vqcr, ~, ~] = solve_opf_qcr(casedata,model);
[rtcr, ctcr, xtcr] = find(Vtcr); [rqcr, cqcr, xqcr] = find(Vqcr);
vtcr = soltcr{1}; vvtcr = vtcr*vtcr'; ytcr = vvtcr(sub2ind([length(vtcr) length(vtcr)],rtcr,ctcr));
vqcr = solqcr{1}; vvqcr = vqcr*vqcr'; yqcr = vvqcr(sub2ind([length(vqcr) length(vqcr)],rqcr,cqcr));
% Exactness error (epsilon)
error = 100*[norm(xtcr-ytcr)/norm(xtcr), norm(xqcr-yqcr)/norm(xqcr)];
[optval, optsol] = solve_opf(casedata,model);
% Optimality gap (gamma)
opt = [opttcr, optqcr];
optgap = 100*abs(1 - opt/max(max(opt),optval));
% Optimality distance (rho)
solgap = 100*[norm(optsol{1}-vtcr)/norm(optsol{1}), norm(optsol{1}-vqcr)/norm(optsol{1})];
isexact = (error < 0.01);