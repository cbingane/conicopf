% 'solve_opf_qcr.m' solves the quadratic convex relaxation of the ACOPF
% problem
function [optval, optsol, Vopt, cpu, status] = solve_opf_qcr(casedata,model)
[n, slack, angslack, pL, qL, gs, bs, vl, vu,...
    nGen, pGl, pGu, qGl, qGu, c2, c1, c0, busgen,...
    nBranch, from, to, y, bsh, tap, shift, su, dl, du,...
    incidentF, incidentT, edges] = opf_data(casedata,model);
Adj = adjacency(graph(edges(:,1),edges(:,2)));
wl = vl(from).*vl(to); wu = vu(from).*vu(to);
csl = -ones(nBranch,1); csu = ones(nBranch,1); snl = csl; snu = csu;
csl(find(dl > -pi/2 & du < pi/2)) = cos(dl(find(dl > -pi/2 & du < pi/2)));
snl(find(dl > -pi/2 & du < pi/2)) = sin(dl(find(dl > -pi/2 & du < pi/2)));
snu(find(dl > -pi/2 & du < pi/2)) = sin(du(find(dl > -pi/2 & du < pi/2)));
Yft = makeYft(nBranch,y,bsh,tap,shift);
cvx_begin
%     cvx_precision low
%     cvx_solver gurobi
    variable V(n,n) hermitian
    variables pf(nBranch) qf(nBranch) pt(nBranch) qt(nBranch)
    variables pG(nGen) qG(nGen)
    variables vm(n) va(n)
    variables w(nBranch) cs(nBranch) sn(nBranch)
    minimize ( c2'*(pG.^2) + c1'*pG + sum(c0) )
    subject to
        for k=1:n
            % PF EQUATIONS
            busgen(:,k)'*pG - pL(k) - gs(k)*V(k,k) == incidentF(:,k)'*pf + incidentT(:,k)'*pt
            busgen(:,k)'*qG - qL(k) + bs(k)*V(k,k) == incidentF(:,k)'*qf + incidentT(:,k)'*qt
            % VOLTAGE LIMITS
            vm(k)^2 <= V(k,k) <= (vl(k) + vu(k))*vm(k) - vl(k)*vu(k)
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
            % FLOW LIMITS
            if (su(l) ~= 0)
                abs(pf(l) + 1j*qf(l)) <= su(l)
                abs(pt(l) + 1j*qt(l)) <= su(l)
            end
            % DIFF PHASE LIMITS
            if (dl(l) > -pi/2 && du(l) < pi/2)
                tan(dl(l))*real(V(from(l),to(l))) <= imag(V(from(l),to(l))) <= tan(du(l))*real(V(from(l),to(l)))
                % cs_l = cos(va_k - va_m), sn_l = sin(va_k - va_m)
                dl(l) <= va(from(l)) - va(to(l)) <= du(l)
                cs(l) >= csl(l)
                cs(l) <= 1 - ((1 - csl(l))/du(l)^2)*(va(from(l)) - va(to(l)))^2
                sn(l) <= cos(du(l)/2)*((va(from(l)) - va(to(l))) - du(l)/2) + sin(du(l)/2)
                sn(l) >= cos(du(l)/2)*((va(from(l)) - va(to(l))) + du(l)/2) - sin(du(l)/2)
                % tightening
                real(V(from(l),to(l)))*cos((du(l)+dl(l))/2) + imag(V(from(l),to(l)))*sin((du(l)+dl(l))/2) >= (vl(from(l))*vl(to(l)) + (vl(to(l))/(vl(from(l)) + vu(from(l))))*(V(from(l),from(l)) - vl(from(l))^2) + (vl(from(l))/(vl(to(l)) + vu(to(l))))*(V(to(l),to(l)) - vl(to(l))^2))*cos((du(l)-dl(l))/2)
                real(V(from(l),to(l)))*cos((du(l)+dl(l))/2) + imag(V(from(l),to(l)))*sin((du(l)+dl(l))/2) >= (vu(from(l))*vu(to(l)) - (vu(to(l))/(vl(from(l)) + vu(from(l))))*(vu(from(l))^2 - V(from(l),from(l))) - (vu(from(l))/(vl(to(l)) + vu(to(l))))*(vu(to(l))^2 - V(to(l),to(l))))*cos((du(l)-dl(l))/2)
            end
            % QCR
            % w_l = v_k*v_m
            w(l) >= vm(from(l))*vl(to(l)) + vl(from(l))*vm(to(l)) - vl(from(l))*vl(to(l))
            w(l) >= vm(from(l))*vu(to(l)) + vu(from(l))*vm(to(l)) - vu(from(l))*vu(to(l))
            w(l) <= vm(from(l))*vu(to(l)) + vl(from(l))*vm(to(l)) - vl(from(l))*vu(to(l))
            w(l) <= vm(from(l))*vl(to(l)) + vu(from(l))*vm(to(l)) - vu(from(l))*vl(to(l))
            % Re(V_km) = w_l*cs_l
            real(V(from(l),to(l))) >= w(l)*csl(l) + wl(l)*cs(l) - wl(l)*csl(l)
            real(V(from(l),to(l))) >= w(l)*csu(l) + wu(l)*cs(l) - wu(l)*csu(l)
            real(V(from(l),to(l))) <= w(l)*csu(l) + wl(l)*cs(l) - wl(l)*csu(l)
            real(V(from(l),to(l))) <= w(l)*csl(l) + wu(l)*cs(l) - wu(l)*csl(l)
            % Im(V_km) = w_l*sn_l
            imag(V(from(l),to(l))) >= w(l)*snl(l) + wl(l)*sn(l) - wl(l)*snl(l)
            imag(V(from(l),to(l))) >= w(l)*snu(l) + wu(l)*sn(l) - wu(l)*snu(l)
            imag(V(from(l),to(l))) <= w(l)*snu(l) + wl(l)*sn(l) - wl(l)*snu(l)
            imag(V(from(l),to(l))) <= w(l)*snl(l) + wu(l)*sn(l) - wu(l)*snl(l)
            % cs_l^2 + sn_l^2 = 1
            abs(cs(l) + 1j*sn(l)) <= 1
        end
        % SLACK BUS
        va(slack) == angslack
cvx_end
% Optimal solution
optval = cvx_optval; Vopt = V; cpu = cvx_cputime; status = cvx_status;
optsol = {approx_volt_profile(Adj,Vopt,slack,angslack); pG + 1j*qG};
end