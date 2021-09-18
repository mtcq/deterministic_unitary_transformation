function P=basis_Uk_overcomplete_method(d,k)
%Creates a basis for the commutant of U^{\otimes k} based on permutation
%operators

%Start by creating a set of k! permutation operators
list_of_perm=perms([1:k]);
V=nan(d^(k),d^(k),factorial(k));
for i=1:factorial(k)
    aux=list_of_perm(i,:);
   V(:,:,i)=PermutationOperator(d,aux);
end

%Write all these operators as vectors to use the GramSchmidit algorithm
for i=1:factorial(k)
   aux=V(:,:,i);
   aux=aux(:);
   V_vectors(:,i)=aux;
end

%Peform the GramSchmidit algorithm
rank_vec=rank(V_vectors);
V_vectors_GRAM=GramSchmidtOrtonomalisation(V_vectors);

P=nan(d^(k),d^(k),rank_vec);
for i=1:rank_vec
    aux=V_vectors_GRAM(:,i);
    aux=reshape(aux,d^(k),d^(k));
   P(:,:,i) = aux;
end

end