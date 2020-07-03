function [optval,optsol,time] = solve_orpd(casedata,model,u,t)
casedata = loadcase(casedata);
mpc = ext2int(casedata);
if (isempty(u) && isempty(t))
    [optval,optsol,time] = solve_opf(mpc,model);
else
    if (~isempty(u) && isempty(t))
        u(:,2) = round(u(:,2));
        mpc.bus(u(:,1),5) = mpc.bus(u(:,1),5).*u(:,2);
        mpc.bus(u(:,1),6) = mpc.bus(u(:,1),6).*u(:,2);
        [optval,optsol,time] = solve_opf(mpc,model);
    end
    if (isempty(u) && ~isempty(t))
        t(:,2) = 1 + round((t(:,2) - 1)/0.0125)*0.0125;
        mpc.branch(t(:,1),9) = t(:,2);
        [optval,optsol,time] = solve_opf(mpc,model);
    end
    if (~isempty(u) && ~ isempty(t))
        u(:,2) = round(u(:,2));
        t(:,2) = 1 + round((t(:,2) - 1)/0.0125)*0.0125;
        mpc.bus(u(:,1),5) = mpc.bus(u(:,1),5).*u(:,2);
        mpc.bus(u(:,1),6) = mpc.bus(u(:,1),6).*u(:,2);
        mpc.branch(t(:,1),9) = t(:,2);
        [optval,optsol,time] = solve_opf(mpc,model);
    end
end
end