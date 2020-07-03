% 'solve_mpopf_qcr.m' solves the quadratic convex relaxation of the
% multi-period ACOPF problem
% INPUTS
%   casedata: MATPOWER case
%   model: either 0 for loss minimization or 1 for cost minimization
function [optval, optsol, Vopt, cpu, status] = solve_mpopf_qcr(casedata,model)
[n, slack, angslack, pL, qL, gs, bs, vl, vu,...
    nGen, pGl, pGu, qGl, qGu, c2, c1, c0, busgen,...
    nBranch, from, to, y, bsh, tap, shift, su, dl, du,...
    incidentF, incidentT, edges] = opf_data(casedata, model);
Adj = adjacency(graph(edges(:,1),edges(:,2)));
wl = vl(from).*vl(to); wu = vu(from).*vu(to);
csl = -ones(nBranch,1); csu = ones(nBranch,1); snl = csl; snu = csu;
csl(find(dl > -pi/2 & du < pi/2)) = cos(dl(find(dl > -pi/2 & du < pi/2)));
snl(find(dl > -pi/2 & du < pi/2)) = sin(dl(find(dl > -pi/2 & du < pi/2)));
snu(find(dl > -pi/2 & du < pi/2)) = sin(du(find(dl > -pi/2 & du < pi/2)));
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
    variable v(n,Tp) complex
    variables pf(nBranch,Tp) qf(nBranch,Tp) pt(nBranch,Tp) qt(nBranch,Tp)
    variables pG(nGen,Tp) qG(nGen,Tp)
    variables vm(n,Tp) va(n,Tp)
    variables w(nBranch,Tp) cs(nBranch,Tp) sn(nBranch,Tp)
    minimize ( sum(c2'*(pG.^2) + c1'*pG + sum(c0)) )
    subject to
        for tp=1:Tp
            for k=1:n
                % PF EQUATIONS
                busgen(:,k)'*pG(:,tp) - pL(k,tp) - gs(k)*V(k,k,tp) == incidentF(:,k)'*pf(:,tp) + incidentT(:,k)'*pt(:,tp)
                busgen(:,k)'*qG(:,tp) - qL(k,tp) + bs(k)*V(k,k,tp) == incidentF(:,k)'*qf(:,tp) + incidentT(:,k)'*qt(:,tp)
                % VOLTAGE LIMITS
                vm(k,tp)^2 <= V(k,k,tp) <= (vl(k) + vu(k))*vm(k,tp) - vl(k)*vu(k)
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
                % SOCR
                [V(from(l),from(l),tp) V(from(l),to(l),tp);...
                    V(to(l),from(l),tp) V(to(l),to(l),tp)] == hermitian_semidefinite(2)
                % FLOW LIMITS
                if (su(l) ~= 0)
                    pf(l,tp)^2 + qf(l,tp)^2 <= su(l)^2
                    pt(l,tp)^2 + qt(l,tp)^2 <= su(l)^2
                end
                % DIFF PHASE LIMITS
                if (dl(l) > -pi/2 && du(l) < pi/2)
                    tan(dl(l))*real(V(from(l),to(l),tp)) <= imag(V(from(l),to(l),tp)) <= tan(du(l))*real(V(from(l),to(l),tp))
                    % cs_l = cos(d_k - d_m), sn_l = sin(d_k - d_m)
                    dl(l) <= va(from(l),tp) - va(to(l),tp) <= du(l)
                    cs(l,tp) >= csl(l)
                    cs(l,tp) <= 1 - ((1 - csl(l))/du(l)^2)*(va(from(l),tp) - va(to(l),tp))^2
                    sn(l,tp) <= cos(du(l)/2)*((va(from(l),tp) - va(to(l),tp)) - du(l)/2) + sin(du(l)/2)
                    sn(l,tp) >= cos(du(l)/2)*((va(from(l),tp) - va(to(l),tp)) + du(l)/2) - sin(du(l)/2)
                end
                % QCR
                % w_l = v_k*v_m
                w(l,tp) >= vm(from(l),tp)*vl(to(l)) + vl(from(l))*vm(to(l),tp) - vl(from(l))*vl(to(l))
                w(l,tp) >= vm(from(l),tp)*vu(to(l)) + vu(from(l))*vm(to(l),tp) - vu(from(l))*vu(to(l))
                w(l,tp) <= vm(from(l),tp)*vu(to(l)) + vl(from(l))*vm(to(l),tp) - vl(from(l))*vu(to(l))
                w(l,tp) <= vm(from(l),tp)*vl(to(l)) + vu(from(l))*vm(to(l),tp) - vu(from(l))*vl(to(l))
                % Re(V_km) = w_l*cs_l
                real(V(from(l),to(l),tp)) >= w(l,tp)*csl(l) + wl(l)*cs(l,tp) - wl(l)*csl(l)
                real(V(from(l),to(l),tp)) >= w(l,tp)*csu(l) + wu(l)*cs(l,tp) - wu(l)*csu(l)
                real(V(from(l),to(l),tp)) <= w(l,tp)*csu(l) + wl(l)*cs(l,tp) - wl(l)*csu(l)
                real(V(from(l),to(l),tp)) <= w(l,tp)*csl(l) + wu(l)*cs(l,tp) - wu(l)*csl(l)
                % Im(V_km) = w_l*sn_l
                imag(V(from(l),to(l),tp)) >= w(l,tp)*snl(l) + wl(l)*sn(l,tp) - wl(l)*snl(l)
                imag(V(from(l),to(l),tp)) >= w(l,tp)*snu(l) + wu(l)*sn(l,tp) - wu(l)*snu(l)
                imag(V(from(l),to(l),tp)) <= w(l,tp)*snu(l) + wl(l)*sn(l,tp) - wl(l)*snu(l)
                imag(V(from(l),to(l),tp)) <= w(l,tp)*snl(l) + wu(l)*sn(l,tp) - wu(l)*snl(l)
                % cs_l^2 + sn_l^2 = 1
                cs(l,tp)^2 + sn(l,tp)^2 <= 1
            end
            % SLACK BUS
            va(slack,tp) == angslack
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