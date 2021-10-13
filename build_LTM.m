function LTM = build_LTM(m,Hp)
% Constructs lower triangular matrix of identity matrices of size m x m.

I = eye(m);
LTM = zeros(m*Hp, m*Hp);

for i = 1:Hp
    LTM((i-1)*m+1:i*m, 1:m*i) = repmat(I,1,i);
end