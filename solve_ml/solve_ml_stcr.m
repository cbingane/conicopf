% 'solve_ml_stcr.m' solves the strong tight-and-cheap relaxation of the ML
% problem
function [optval, optsol, Vopt, cpu, status] = solve_ml_stcr(casedata)
[n, slack, angslack, pL, qL, gs, bs, vl, vu,...
    nGen, pGl, pGu, qGl, qGu, ~, ~, ~, busgen,...
    nBranch, from, to, y, bsh, tap, shift, su, dl, du,...
    incidentF, incidentT, edges] = opf_data(casedata);
Adj = adjacency(graph(edges(:,1),edges(:,2)));
Yft = makeYft(nBranch,y,bsh,tap,shift);
cvx_begin
%     cvx_precision low
%     cvx_solver sdpt3
    variable V(n,n) hermitian
    variables pf(nBranch) qf(nBranch) pt(nBranch) qt(nBranch)
    variables pG(nGen) qG(nGen)
    variable ml
    maximize ( ml )
    subject to
        ml >= 1
        for k=1:n
            % PF EQUATIONS
            busgen(:,k)'*pG - pL(k)*ml - gs(k)*V(k,k) == incidentF(:,k)'*pf + incidentT(:,k)'*pt
            busgen(:,k)'*qG - qL(k)*ml + bs(k)*V(k,k) == incidentF(:,k)'*qf + incidentT(:,k)'*qt
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
            % STCR
            if (~ismember(slack,[from(l) to(l)]))
                [V(slack,slack) V(slack,from(l)) V(slack,to(l));...
                    V(from(l),slack) V(from(l),from(l)) V(from(l),to(l));...
                    V(to(l),slack) V(to(l),from(l)) V(to(l),to(l))] == hermitian_semidefinite(3)
            elseif ((to(l) == slack && degree(graph(Adj),from(l)) == 1) || (from(l) == slack && degree(graph(Adj),to(l)) == 1))
               [V(from(l),from(l)) V(from(l),to(l));...
                   V(to(l),from(l)) V(to(l),to(l))] == hermitian_semidefinite(2) 
            end
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