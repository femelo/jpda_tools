function I = intersect2(A,B)

na = numel(A); % A = [1 2 3]; na = 3
nb = numel(B); % B = [1 4]; nb = 2

if na > nb
    Ca = ones(nb,1) * A;    % Ca = [1 2 3; 1 2 3];
    Cb = (B') * ones(1,na); % Cb = [1 1 1; 4 4 4];
    D = prod(Ca-Cb,2)';
    I = B(D == 0);
else
    Cb = ones(na,1) * B;    % Ca = [1 2 3; 1 2 3];
    Ca = (A') * ones(1,nb); % Cb = [1 1 1; 4 4 4];
    D = prod(Ca-Cb,2)';
    I = A(D == 0);
end

return;