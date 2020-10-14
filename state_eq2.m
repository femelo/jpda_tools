function [xk, Pk] = state_eq2(xkm1, Pkm1, ukm1, Qkm1)

T = 1;

F2 = [1 T; 0 1];
% G2 = [T^2/2 0; 0 T];
G2 = eye(2);

Fkm1 = blkdiag(F2,F2);
Gkm1 = blkdiag(G2,G2);

xk = Fkm1*xkm1 + Gkm1*ukm1;
Pk = Fkm1*Pkm1*Fkm1' + Qkm1;

end

