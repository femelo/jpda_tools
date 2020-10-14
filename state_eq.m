function [xk, Pk] = state_eq(xkm1, Pkm1, ukm1, Qkm1)

Fkm1 = eye(2);

xk = Fkm1*xkm1 + ukm1;
Pk = Fkm1*Pkm1*Fkm1' + Qkm1;

end

