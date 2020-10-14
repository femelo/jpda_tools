function C = union2(a,B)

D = B-a;
flag = prod(D) == 0;

if flag
   C = B;
else
   C = sort([a B]);
end

return;