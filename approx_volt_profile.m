function v = approx_volt_profile(Adj,V,slack,angslack)
% Approximate voltage profile
n = size(Adj,1);
modv = sqrt(diag(V)); angv = sparse(n,1);
angv(slack) = angslack;
for k = 1:n
    if k ~= slack
        path = shortestpath(graph(Adj),slack,k);
        p = length(path);
        while p > 1
            angv(k) = angv(k) + angle(V(path(p), path(p-1)));
            p = p - 1;
        end
        angv(k) = angv(k) + angslack;
    end
end
v = full(modv.*exp(1j*angv));
end