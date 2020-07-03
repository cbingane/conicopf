function [optval, optsol, Vopt, topt, cpu, status] = solve_orpd_tap_sdr1(casedata,model)
[n, slack, angslack, pL, qL, gs, bs, vl, vu,...
    nGen, pGl, pGu, qGl, qGu, c2, c1, c0, busgen,...
    nBranch, from, to, y, bsh, tap, shift, su, dl, du, incidentF, incidentT, edges] = opf_data(casedata, model);
isTap = find(spones(tap).*(1 - spones(shift)));
if (isempty(isTap))
    [optval, optsol, Vopt, cpu, status] = solve_opf_sdr(casedata,model);
    topt = [];
else
Adj = adjacency(graph(edges(:,1),edges(:,2)));
Yft = makeYft(nBranch,y,bsh,tap,shift);
tl = 0.9*ones(length(isTap),1); tu = 1.1*ones(length(isTap),1);
cvx_begin
%     cvx_precision low
%     cvx_solver sdpt3
    variable V(n,n) hermitian
    variables pf(nBranch) qf(nBranch) pt(nBranch) qt(nBranch)
    variables pG(nGen) qG(nGen)
    variable W(length(isTap))
    variable Wf(length(isTap))
    variable Wt(length(isTap)) complex
    minimize ( c2'*(pG.^2) + c1'*pG + sum(c0) )
    subject to
        % SDR
        V == hermitian_semidefinite(n)
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
            if (ismember(l,isTap))
                % WITH TAP
                pf(l) + 1j*qf(l) == conj(y(l) + 1j*bsh(l)/2)*W(find(isTap == l)) + conj(-y(l))*Wt(find(isTap == l))
                pt(l) + 1j*qt(l) == conj(-y(l))*conj(Wt(find(isTap == l))) + conj(y(l) + 1j*bsh(l)/2)*V(to(l),to(l))
                % CVX FRACTIONAL
                tl(find(isTap == l))^2*W(find(isTap == l)) <= V(from(l),from(l)) <= tu(find(isTap == l))*Wf(find(isTap == l))
                % SDP
                [V(from(l),from(l)) Wf(find(isTap == l)) V(from(l),to(l));...
                    Wf(find(isTap == l)) W(find(isTap == l)) Wt(find(isTap == l));...
                    V(to(l),from(l)) Wt(find(isTap == l))' V(to(l),to(l))] == hermitian_semidefinite(3)
            else
                % NO TAP
                pf(l) + 1j*qf(l) == conj(Yft{l}(1,1))*V(from(l),from(l)) + conj(Yft{l}(1,2))*V(from(l),to(l))
                pt(l) + 1j*qt(l) == conj(Yft{l}(2,1))*V(to(l),from(l)) + conj(Yft{l}(2,2))*V(to(l),to(l))
            end
            % FLOW LIMITS
            if (su(l) ~= 0)
                pf(l)^2 + qf(l)^2 <= su(l)^2
                pt(l)^2 + qt(l)^2 <= su(l)^2
            end
            % DIFF PHASE LIMITS
            if (dl(l) > -pi/2 && du(l) < pi/2)
                tan(dl(l))*real(V(from(l),to(l))) <= imag(V(from(l),to(l))) <= tan(du(l))*real(V(from(l),to(l)))
            end
        end
cvx_end
% Optimal solution
optval = cvx_optval; Vopt = V; cpu = cvx_cputime; status = cvx_status;
optsol = {approx_volt_profile(Adj,Vopt,slack,angslack); pG + 1j*qG};
tW = ones(length(isTap),1);
for i=1:length(isTap)
    tW(i) = sqrt(V(from(isTap(i)),from(isTap(i)))/W(i));
end
topt = [isTap tW];
end
end