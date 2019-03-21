% 'is_tcr_exact.m' measures the exactness of the tight-and-cheap relaxation
% for the ACOPF problem
function [optval, opttcr, optgap, solgap, metric, isexact] = is_tcr_exact(casedata,model)
[Vtcr, vtcr, ~, ~, opttcr] = solve_opf_tcr(casedata,model);
vtcr = vtcr(:,1);
[optval, optsol] = solve_opf(casedata,model);
optgap = 100*(1 - opttcr/max(opttcr,optval));
solgap = 100*norm(optsol - vtcr)/norm(optsol);
% metric = 100*norm(diag(Vtcr) - diag(vtcr*vtcr'))/norm(diag(Vtcr));
error = 100*(1 -  abs(vtcr)./sqrt(diag(Vtcr)));
metric = [min(error) mean(error) max(error)];
isexact = (metric < 0.01);