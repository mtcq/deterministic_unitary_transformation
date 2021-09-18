function A = GramSchmidtOrtonomalisation(A)
%Inputs: a matrix with a set of column vectors (which are not necessarily
%all linearly independent
%Output: An orthonormal basis for the space spanned by vectors in A

%We start by removing all the linearly dependend column vector in A such
%that A is a full rank matrix
r=rank(A);
[Q, R, E] = qr(A,0); 
indices=sort(E(1:r));
A=A(:,indices);

%Perform the standard Gram-Schmidit algorithm
for i = 1:r
    A(:,i) = A(:,i) / norm(A(:,i));
    for j = i+1:r
       A(:,j) = A(:,j) - A(:,i)'*A(:,j) /norm(A(:,i)) * A(:,i);
    end
end

end