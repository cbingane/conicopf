% 'solve_opf_chr.m' solves the chordal relaxation of the ACOPF problem
function [Vchr, vchr, cpuchr, statuschr, optchr] = solve_opf_chr(casedata,model)
[n, slack, angslack, pL, qL, gs, bs, vl, vu,...
    nGen, pGl, pGu, qGl, qGu, c2, c1, c0, busgen,...
    nBranch, from, to, y, bsh, tap, shift, su, dl, du, incidentF, incidentT, edges] = opf_data(casedata, model);
Adj = adjacency(graph(edges(:,1),edges(:,2)));
Vl = vl.^2; Vu = vu.^2;
[~, bag] = chordamd(Adj);
Yft = makeYft(nBranch,y,bsh,tap,shift);
cvx_begin
%     cvx_precision low
%     cvx_solver sedumi
    variable V(n,n) hermitian
    variables pf(nBranch) qf(nBranch) pt(nBranch) qt(nBranch)
    variables pG(nGen) qG(nGen)
    minimize ( c2'*(pG.^2) + c1'*pG + sum(c0) )
    subject to
        % CHR
        for k = 1:n
            if (~isempty(bag{k}))
                V(bag{k},bag{k}) == hermitian_semidefinite(length(bag{k}))
            end
        end
        for k=1:n
            % PF EQUATIONS
            busgen(:,k)'*pG - pL(k) - gs(k)*V(k,k) == incidentF(:,k)'*pf + incidentT(:,k)'*pt
            busgen(:,k)'*qG - qL(k) + bs(k)*V(k,k) == incidentF(:,k)'*qf + incidentT(:,k)'*qt
            % VOLTAGE LIMITS
            Vl(k) <= V(k,k) <= Vu(k)
        end
        % GENERATION LIMITS
        for g = 1:nGen
            pGl(g) <= pG(g) <= pGu(g)
            qGl(g) <= qG(g) <= qGu(g)
        end
        % BRANCH CONSTRAINTS
%         V((Ach + speye(n)) == 0) == 0
        for l=1:nBranch
            % FLOW INJECTION
            pf(l) + 1j*qf(l) == conj(Yft{l}(1,1))*V(from(l),from(l)) + conj(Yft{l}(1,2))*V(from(l),to(l))
            pt(l) + 1j*qt(l) == conj(Yft{l}(2,1))*V(to(l),from(l)) + conj(Yft{l}(2,2))*V(to(l),to(l))
            % FLOW LIMITS
            if (su(l) ~= 0)
                [su(l) pf(l) + 1j*qf(l); pf(l) - 1j*qf(l) su(l)] == hermitian_semidefinite(2)
                [su(l) pt(l) + 1j*qt(l); pt(l) - 1j*qt(l) su(l)] == hermitian_semidefinite(2)
            end
            % DIFF PHASE LIMITS
            if (dl(l) > -pi/2 && du(l) < pi/2)
                tan(dl(l))*real(V(from(l),to(l))) <= imag(V(from(l),to(l))) <= tan(du(l))*real(V(from(l),to(l)))
            end
        end
cvx_end
% Optimal solution
cpuchr = cvx_cputime;
optchr = cvx_optval; statuschr = cvx_status;
Vchr = V;
vchr = approx_volt_profile(Adj,Vchr,slack,angslack);
end