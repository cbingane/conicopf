% 'is_tcr_exact.m' measures the exactness of the tight-and-cheap relaxation
% for the ACOPF problem
function [optval, opttcr, optgap, solgap, error, isexact] = is_tcr_exact(casedata,model)
[Vtcr, vtcr, ~, ~, opttcr] = solve_opf_tcr(casedata,model);
vtcr = vtcr(:,1);
[optval, optsol] = solve_opf(casedata,model);
optgap = 100*(1 - opttcr/max(opttcr,optval));
solgap = 100*norm(optsol - vtcr)/norm(optsol);
error = 100*max(1 -  abs(vtcr)./sqrt(diag(Vtcr)));
isexact = (error < 0.01);