## Code to accompany: *[Deterministic transformations between unitary operations: Exponential advantage with adaptive quantum circuits and the power of indefinite causality](https://arxiv.org/abs/xxx)*
#### Marco Túlio Quintino and Daniel Ebler



 MATLAB code requires:
- [cvx](http://cvxr.com/) - a free MATLAB toolbox for rapid prototyping of optimization problems.
- [QETLAB](http://www.qetlab.com/) - a free MATLAB toolbox for quantum entanglement theory.

 MATHEMATICA code requires:
- [QI](https://github.com/rogercolbeck/QI) - a free quantum information Mathematica package.

The MATLAB code of this repository is:

- [run_max_F.m](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/run_max_F.m):
Script used to evaluate the maximal average fidelity for transforming k copies of a d-dimensional unitary transformations (transposition, complex conjugation and inversion) based on the SDP primal.

- [run_max_F_dual.m](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/run_max_F.m):
Script used to evaluate the maximal average fidelity for transforming k copies of a d-dimensional unitary transformations (transposition, complex conjugation and inversion) based on the SDP dual.

- [store_matlab_variables.m](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/store_matlab_variables.m):
Script used to evaluate the maximal average fidelity for transforming k copies of a d-dimensional unitary transformations (transposition, complex conjugation and inversion) for all cases considered in the paper. This script also creates and stores all variables required for the computer assisted proof.

- [optimal_fU.m](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/optimal_fU.m):
Function that evaluates the maximal average fidelity for transforming k copies of a d-dimensional unitary transformations (transposition, complex conjugation and inversion) based on the SDP dual.

- [make_OMEGA_permutation.m](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/make_OMEGA_permutation.m):
Function uses to construct the performance operator Omega for transforming k copies of a d-dimensional unitary transformations (transposition, complex conjugation and inversion). This function constructs Omega based on the permutation operators.

- [make_OMEGA_group.m](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/make_OMEGA_group.m):
Function uses to construct the performance operator Omega for transforming k copies of a d-dimensional unitary transformations (transposition, complex conjugation and inversion). This function constructs Omega based on the explicit basis obtain via group theory methods.

- [basis_UUU.m](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/basis_UUU.m):
Function that construct a basis for the commutant of U⊗U⊗U  based on group theory methods.

- [basis_UstarUstarU.m](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/basis_UstarUstarU.m):
Function that construct a basis for the commutant of U*⊗U*⊗U* based on group theory methods.

- [basis_Uk_overcomplete_method.m](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/basis_Uk_overcomplete_method.m):
Function that construct a basis for the commutant of k copies of a d-dimension unitary U based on permutation operators.

- [basis_Ukstar_overcomplete_method.m](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/basis_Ukstar_overcomplete_method.m):
Function that construct a basis for the commutant of k copies of a d-dimension unitary U where some are conjugated based on permutation operators.

- [is_general_superchannel.m](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/TR.m):
Function that imposes the SDP constraint of a k-slot general superchannel.

- [GramSchmidtOrtonomalisation.m](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/GramSchmidtOrtonomalisation.m):
Function that performs the Gram-Schmidt orthonormalisation algorithm on a set of vectors.

- [TR.m](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/TR.m):
Function that implements the trace-and-replace map on matrices.

The MATHEMATICA code of this repository is:

- [LoadFunctions.nb](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/LoadFunctions.nb):
Notebook used to load important functions.

- [ComputerAssitedProof.nb](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/ComputerAssitedProof.nb):
Notebook used to certify the qubit values for unitary transposition in the table of the paper.

- [ComputerAssitedProofd3.nb](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/ComputerAssitedProofd3.nb):
Notebook used to certify the qutrit values for unitary transposition, presented in the table of the paper.

- [ComputerAssitedProofd3INV.nb](https://github.com/mtcq/deterministic_unitary_transformation/blob/main/code/ComputerAssitedProofd3INV.nb):
Notebook used to certify the qutrit values for unitary inversion presented in the table of the paper.
