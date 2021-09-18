function P=basis_UstarUstarU(d)

dketbraId=IsotropicState(d,1)*d;
F=SwapOperator(d);
%X=PermuteSystems(kron(dketbraId,eye(d)),[3 2 1],[d d d]);
X=kron(eye(d),dketbraId);
%X=kron(dketbraId,eye(d));
%V=PermuteSystems(kron(eye(d),F),[3 2 1],[d d d]);
V=kron(F,eye(d));
Id123=eye(d^3);

Sp = (Id123+V)/2*(Id123-2*X/(d+1))*(Id123+V)/2;
Sm = (Id123-V)/2*(Id123-2*X/(d-1))*(Id123-V)/2;
S0 = 1/(d^2-1)*(d*(X+V*X*V)-(X*V+V*X));

S1 = 1/(d^2-1)*(d*(X*V+V*X)-(X+V*X*V));
S2 = 1/sqrt(d^2-1)*(X-V*X*V);
S3 = sqrt(-1)/sqrt(d^2-1)*(X*V-V*X);


P(:,:,1) = Sp;
P(:,:,2) = Sm;
P(:,:,3) = S0;

P(:,:,4) = S1;
P(:,:,5) = S2;
P(:,:,6) = S3;

end
