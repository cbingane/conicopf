% 'solve_opf_socr.m' solves the second-order cone relaxation of the ACOPF
% problem
function [optval, optsol, Vopt, cpu, status] = solve_opf_socr(casedata,model)
[n, slack, angslack, pL, qL, gs, bs, vl, vu,...
    nGen, pGl, pGu, qGl, qGu, c2, c1, c0, busgen,...
    nBranch, from, to, y, bsh, tap, shift, su, dl, du,...
    incidentF, incidentT, edges] = opf_data(casedata,model);
Adj = adjacency(graph(edges(:,1),edges(:,2)));
Yft = makeYft(nBranch,y,bsh,tap,shift);
cvx_begin
%     cvx_precision low
%     cvx_solver gurobi
    variable V(n,n) hermitian
    variables pf(nBranch) qf(nBranch) pt(nBranch) qt(nBranch)
    variables pG(nGen) qG(nGen)
    minimize ( c2'*(pG.^2) + c1'*pG + sum(c0) )
    subject to
        for k=1:n
            % PF EQUATIONS
            busgen(:,k)'*pG - pL(k) - gs(k)*V(k,k) == incidentF(:,k)'*pf + incidentT(:,k)'*pt
            busgen(:,k)'*qG - qL(k) + bs(k)*V(k,k) == incidentF(:,k)'*qf + incidentT(:,k)'*qt
            % VOLTAGE LIMITS
            vl(k)^2 <= V(k,k) <= vu(k)^2
        end
        % GENERATION LIMITS
        for g = 1:nGen
            pGl(g) <= pG(g) <= pGu(g)
            qGl(g) <= qG(g) <= qGu(g)
        end
        % BRANCH CONSTRAINTS
        for l=1:nBranch
            % FLOW INJECTION
            pf(l) + 1j*qf(l) == conj(Yft{l}(1,1))*V(from(l),from(l)) + conj(Yft{l}(1,2))*V(from(l),to(l))
            pt(l) + 1j*qt(l) == conj(Yft{l}(2,1))*V(to(l),from(l)) + conj(Yft{l}(2,2))*V(to(l),to(l))
            % SOCR
            [V(from(l),from(l)) V(from(l),to(l)); V(to(l),from(l)) V(to(l),to(l))] == hermitian_semidefinite(2)
%             [V(from(l),from(l)) + V(to(l),to(l)) 0 2*V(from(l),to(l));...
%                 0 V(from(l),from(l)) + V(to(l),to(l)) V(from(l),from(l)) - V(to(l),to(l));...
%                 2*V(to(l),from(l)) V(from(l),from(l)) - V(to(l),to(l)) V(from(l),from(l)) + V(to(l),to(l))] == hermitian_semidefinite(3)
            % FLOW LIMITS
            if (su(l) ~= 0)
                abs(pf(l) + 1j*qf(l)) <= su(l)
                abs(pt(l) + 1j*qt(l)) <= su(l)
            end
            % DIFF PHASE LIMITS
            if (dl(l) > -pi/2 && du(l) < pi/2)
                tan(dl(l))*real(V(from(l),to(l))) <= imag(V(from(l),to(l))) <= tan(du(l))*real(V(from(l),to(l)))
                % tightening
                real(V(from(l),to(l)))*cos((du(l)+dl(l))/2) + imag(V(from(l),to(l)))*sin((du(l)+dl(l))/2) >= (vl(from(l))*vl(to(l)) + (vl(to(l))/(vl(from(l)) + vu(from(l))))*(V(from(l),from(l)) - vl(from(l))^2) + (vl(from(l))/(vl(to(l)) + vu(to(l))))*(V(to(l),to(l)) - vl(to(l))^2))*cos((du(l)-dl(l))/2)
                real(V(from(l),to(l)))*cos((du(l)+dl(l))/2) + imag(V(from(l),to(l)))*sin((du(l)+dl(l))/2) >= (vu(from(l))*vu(to(l)) - (vu(to(l))/(vl(from(l)) + vu(from(l))))*(vu(from(l))^2 - V(from(l),from(l))) - (vu(from(l))/(vl(to(l)) + vu(to(l))))*(vu(to(l))^2 - V(to(l),to(l))))*cos((du(l)-dl(l))/2)
            end
        end
cvx_end
% Optimal solution
optval = cvx_optval; Vopt = V; cpu = cvx_cputime; status = cvx_status;
optsol = {approx_volt_profile(Adj,Vopt,slack,angslack); pG + 1j*qG};
end