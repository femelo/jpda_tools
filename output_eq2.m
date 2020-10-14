function [zk, Sk, Wk] = output_eq2(xk, Pk, vk, Rk)

Hk = [1 0 0 0; 0 0 1 0];

zk = Hk*xk + vk;
Sk = Hk*Pk*Hk' + Rk;

if Sk == 0
    Wk = 0;
else
    Wk = Pk*Hk'*pinv(Sk,eps);
end

end

