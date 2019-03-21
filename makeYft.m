% 'makeYft.m' computes for each branch the elements of the branch
% admittance matrix
function Yft = makeYft(nBranch,y,bsh,tap,shift)
t = ones(nBranch,1);
idxtap = find(tap);
t(idxtap) = tap(idxtap);
tap = t.*exp(1j*shift);
Yft = cell(nBranch,1);
for l = 1:nBranch
    Yft{l}(2,2) = y(l) + 1j*bsh(l)/2;
    Yft{l}(1,1) = Yft{l}(2,2)/abs(tap(l))^2;
    Yft{l}(1,2) = -y(l)/conj(tap(l));
    Yft{l}(2,1) = -y(l)/tap(l);
end
end