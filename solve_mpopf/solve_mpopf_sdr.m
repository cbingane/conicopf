% 'solve_mpopf_sdr.m' solves the semidefinite relaxation of the
% multi-period ACOPF problem
% INPUTS
%   casedata: MATPOWER case
%   model: either 0 for loss minimization or 1 for cost minimization
function [optval, optsol, Vopt, cpu, status] = solve_mpopf_sdr(casedata,model)
[n, slack, angslack, pL, qL, gs, bs, vl, vu,...
    nGen, pGl, pGu, qGl, qGu, c2, c1, c0, busgen,...
    nBranch, from, to, y, bsh, tap, shift, su, dl, du,...
    incidentF, incidentT, edges] = opf_data(casedata, model);
Adj = adjacency(graph(edges(:,1),edges(:,2)));
Yft = makeYft(nBranch,y,bsh,tap,shift);
Tp = 24;
lf = [0.950, 0.953, 0.950, 0.956, 0.962, 1.010, 1.100, 1.052,...
    1.010, 0.995, 0.989, 0.980, 0.980, 0.986, 1.034, 1.073,...
    1.085, 1.067, 1.052, 1.010, 0.995, 0.980, 0.965, 0.956];
pL = pL*ones(Tp,1)'*diag(lf); qL = qL*ones(Tp,1)'*diag(lf);
ramp = sparse(Tp-1,1); rl = ramp; ru = ramp;
for tp=1:Tp-1
    ramp(tp) = sum(pL(:,tp+1) - pL(:,tp))/nGen;
    rl(tp) = ramp(tp) - 0.1*abs(ramp(tp));
    ru(tp) = ramp(tp) + 0.1*abs(ramp(tp));
end
cvx_begin
%     cvx_precision low
%     cvx_solver sdpt3
    variable V(n,n,Tp) hermitian
    variables pf(nBranch,Tp) qf(nBranch,Tp) pt(nBranch,Tp) qt(nBranch,Tp)
    variables pG(nGen,Tp) qG(nGen,Tp)
    minimize ( sum(c2'*(pG.^2) + c1'*pG + sum(c0)) )
    subject to
        for tp=1:Tp
            % SDR
            V(:,:,tp) == hermitian_semidefinite(n)
            for k=1:n
                % PF EQUATIONS
                busgen(:,k)'*pG(:,tp) - pL(k,tp) - gs(k)*V(k,k,tp) == incidentF(:,k)'*pf(:,tp) + incidentT(:,k)'*pt(:,tp)
                busgen(:,k)'*qG(:,tp) - qL(k,tp) + bs(k)*V(k,k,tp) == incidentF(:,k)'*qf(:,tp) + incidentT(:,k)'*qt(:,tp)
                % VOLTAGE LIMITS
                vl(k)^2 <= V(k,k,tp) <= vu(k)^2
            end
            % GENERATION LIMITS
            for g = 1:nGen
                pGl(g) <= pG(g,tp) <= pGu(g)
                qGl(g) <= qG(g,tp) <= qGu(g)
            end
            % BRANCH CONSTRAINTS
            for l=1:nBranch
                % FLOW INJECTION
                pf(l,tp) + 1j*qf(l,tp) == conj(Yft{l}(1,1))*V(from(l),from(l),tp) + conj(Yft{l}(1,2))*V(from(l),to(l),tp)
                pt(l,tp) + 1j*qt(l,tp) == conj(Yft{l}(2,1))*V(to(l),from(l),tp) + conj(Yft{l}(2,2))*V(to(l),to(l),tp)
                % FLOW LIMITS
                if (su(l) ~= 0)
                    pf(l,tp)^2 + qf(l,tp)^2 <= su(l)^2
                    pt(l,tp)^2 + qt(l,tp)^2 <= su(l)^2
                end
                % DIFF PHASE LIMITS
                if (dl(l) > -pi/2 && du(l) < pi/2)
                    tan(dl(l))*real(V(from(l),to(l),tp)) <= imag(V(from(l),to(l),tp)) <= tan(du(l))*real(V(from(l),to(l),tp))
                end
            end
        end
        for tp=1:Tp-1
            rl(tp)*ones(nGen,1) <= pG(:,tp+1) - pG(:,tp) <= ru(tp)*ones(nGen,1)
        end
cvx_end
% Optimal solution
optval = cvx_optval; Vopt = V; cpu = cvx_cputime; status = cvx_status;
optsol = cell(Tp,1);
for tp=1:Tp
    optsol{tp} = {approx_volt_profile(Adj,Vopt(:,:,tp),slack,angslack); pG(:,tp) + 1j*qG(:,tp)};
end
end