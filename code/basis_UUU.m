function P=basis_UUU(d)
%Creates a basis for the commutant of U\otimes U\otimes U based on group
%theory methods

R1 = 1/3*(2*PermutationOperator(d,[1 3 2]) - PermutationOperator(d,[3 2 1]) - PermutationOperator(d,[2 1 3]));
R2 = 1/sqrt(3)*(PermutationOperator(d,[2 1 3]) - PermutationOperator(d,[3 2 1]));
R3 = sqrt(-1)/sqrt(3)*(PermutationOperator(d,[3 1 2]) - PermutationOperator(d,[2 3 1]));

Rp = 1/6*(eye(d^3) + PermutationOperator(d,[2 1 3]) + PermutationOperator(d,[1 3 2]) + PermutationOperator(d,[3 2 1]) + PermutationOperator(d,[3 1 2]) + PermutationOperator(d,[2 3 1]));
Rm = 1/6*(eye(d^3) - PermutationOperator(d,[2 1 3]) - PermutationOperator(d,[1 3 2]) - PermutationOperator(d,[3 2 1]) + PermutationOperator(d,[3 1 2]) + PermutationOperator(d,[2 3 1]));
R0 = 1/3*(2*eye(d^3) - PermutationOperator(d,[3 1 2]) - PermutationOperator(d,[2 3 1]));


P(:,:,1) = R1;
P(:,:,2) = R2;
P(:,:,3) = R3;

P(:,:,4) = Rp;
P(:,:,5) = Rm;
P(:,:,6) = R0;

end
