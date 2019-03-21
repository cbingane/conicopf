function [Ach, bag] = chordamd(Adj)
n = size(Adj,1);
Q = laplacian(graph(Adj)) + speye(n);
idx = amd(Q);
Lchol = chol(Q(idx,idx),'lower');
Asortch = adjacency(graph(Lchol - diag(diag(Lchol)),'lower'));
idxinv(idx) = 1:length(idx);
Ach = Asortch(idxinv,idxinv);
bag = cell(n,1);
I = speye(n); B = Ach; k = 1; nod = numnodes(graph(B));
while (numedges(graph(B)) < nchoosek(nod,2))
    bag{idx(k)} = find(B(:,idx(k))+ I(:,idx(k)));
    B(:,idx(k)) = zeros(n,1);
    B(idx(k),:) = B(:,idx(k))';
    nod = nod-1;
    k = k+1;
end
bag{idx(k)} = find(B(:,idx(k))+ I(:,idx(k)));