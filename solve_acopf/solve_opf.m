% 'solve_opf.m' solves the ACOPF problem and provides the optimal value and
% the optimal bus voltages
% INPUTS
%   casedata: MATPOWER case
%   model: either 0 for loss minimization or 1 for cost minimization
function [optval,optsol,time] = solve_opf(casedata,model)
casedata = loadcase(casedata);
mpc = ext2int(casedata);
mpopt = mpoption('opf.ac.solver','mips');
if model == 0
    % LOSS MINIMIZATION
    if size(mpc.gencost,2) == 7
        mpc.gencost(:,5) = 0;
        mpc.gencost(:,6) = 1;
        mpc.gencost(:,7) = 0;
    end
    if size(mpc.gencost,2) == 6
        mpc.gencost(:,5) = 1;
        mpc.gencost(:,6) = 0;
    end
    sol = runopf(mpc,mpopt);
else
    % COST MINIMIZATION
    sol = runopf(mpc,mpopt);
end
optval = sol.f;
optsol = {sol.x(size(mpc.bus,1)+1:2*size(mpc.bus,1)).*exp(1j*sol.x(1:size(mpc.bus,1)));...
    sol.x(2*size(mpc.bus,1)+1:2*size(mpc.bus,1)+size(mpc.gen,1))+...
    1j*sol.x(2*size(mpc.bus,1)+size(mpc.gen,1)+1:2*size(mpc.bus,1)+2*size(mpc.gen,1))};
time = sol.et;
end