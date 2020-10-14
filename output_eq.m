function [zk, Sk, Wk] = output_eq(xk, Pk, vk, Rk)

Hk = eye(2);

zk = Hk*xk + vk;
Sk = Hk*Pk*Hk' + Rk;

if Sk == 0
    Wk = 0;
else
    Wk = Pk*Hk'*pinv(Sk);
end

end

