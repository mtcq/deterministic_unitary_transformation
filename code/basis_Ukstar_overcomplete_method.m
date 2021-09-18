function P=basis_Ukstar_overcomplete_method(d,k,staredSPACES)
%Creates a basis for the commutant of multiple copies of U and conj(U) based on permutation
%operators and transpositions
%indicate the spaces which which have a complex conjugate in the mask.
% For instance, U*\otimes U* \otimes U has mask [1 2]

list_of_perm=perms([1:k]);
V=nan(d^(k),d^(k),factorial(k));
DIMS=d*ones(1,k);
for i=1:factorial(k)
    aux=list_of_perm(i,:);
   PermutationOperator(d,aux);
   size(PermutationOperator(d,aux));
   V(:,:,i)=PartialTranspose(PermutationOperator(d,aux),staredSPACES,DIMS);
end

for i=1:factorial(k)
   aux=V(:,:,i);
   aux=aux(:);
   V_vectors(:,i)=aux;
end

rank_vec=rank(V_vectors);
V_vectors_GRAM=GramSchmidtOrtonomalisation(V_vectors);

P=nan(d^(k),d^(k),rank_vec);
for i=1:rank_vec
    aux=V_vectors_GRAM(:,i);
    aux=reshape(aux,d^(k),d^(k));
   P(:,:,i) = aux;
end

end