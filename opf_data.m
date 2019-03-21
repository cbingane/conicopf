% 'opf_data.m' reads MATPOWER data file format
function [n, slack, angslack, pL, qL, gs, bs, vl, vu,...
    nGen, pGl, pGu, qGl, qGu, c2, c1, c0, busgen,...
    nBranch, from, to, y, bsh, tap, shift, su, dl, du, incidentF, incidentT, edges] = opf_data(casedata, model)
casedata = loadcase(casedata);
mpc = ext2int(casedata);
%% baseMVA
base = mpc.baseMVA;
%% bus data
bus = mpc.bus; n = size(bus,1);
% slack
slack = find(bus(:,2) == 3); angslack = bus(slack,9)*pi/180;
% loads
pL = sparse(bus(:,3)/base); qL = sparse(bus(:,4)/base);
% bounds Voltages
vl = bus(:,13); vu = bus(:,12);
% shunt compensators
gs = sparse(bus(:,5)/base); bs = sparse(bus(:,6)/base);
%
%% branch data
branch = mpc.branch; nBranch = size(branch,1);
% from to
from = branch(:,1); to = branch(:,2);
% series admittance
y = 1./(branch(:,3) + 1j*branch(:,4));
% shunt conductance
bsh = branch(:,5);
% tap
tap = sparse(branch(:,9));
% shift
shift = sparse(branch(:,10)*pi/180);
% flow limit
su = sparse(branch(:,6)/base);
% flow limit
dl = sparse(branch(:,12)*pi/180); du = sparse(branch(:,13)*pi/180);
% incidence
incidentF = sparse(1:nBranch,from,ones(nBranch,1),nBranch,n);
incidentT = sparse(1:nBranch,to,ones(nBranch,1),nBranch,n);
% adjacence
edges = unique(sort([from to],2), 'rows');
%
%% generator bounds
gen = mpc.gen; nGen = size(gen,1);
pGl = gen(:,10)/base; pGu = gen(:,9)/base;
qGl = gen(:,5)/base; qGu = gen(:,4)/base;
busgen = sparse(1:nGen,gen(:,1),ones(nGen,1),nGen,n);
%% generator cost coefficients
gencost = mpc.gencost;
if size(gencost,2) == 7
    c2 = sparse(gencost(:,5))*base^2; c1 = sparse(gencost(:,6))*base; c0 = sparse(gencost(:,7));
end
if size(gencost,2) == 6
    c2 = sparse(nGen,1)*base^2; c1 = sparse(gencost(:,5))*base; c0 = sparse(gencost(:,6));
end
if model == 0
    c2 = sparse(nGen,1)*base^2;
    c1 = ones(nGen,1)*base;
    c0 = sparse(nGen,1);
end