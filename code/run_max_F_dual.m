clear 
tic
%This is a script that evaluates the maximal average fidelity for unitary
%transformations by using the dual method.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Start of Adjustable Parameters  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d=2               %Dimension 'd'
        k=1              %Number of uses 'k'
        protocol=1; %1 for parallel, 2 for sequential, 3 for general
        type=1;      %1 for Transposition, 2 for Conjugate, 3 for Inversion 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   End of Adjustable Parameters   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
OMEGA= make_OMEGA_permutation(d,k,type);
%OMEGA=make_OMEGA_group(d,k);

switch type
    case 1
        type='transposition'
    case 2
        type='conjugate'
    case 3
        type='inversion'
    otherwise
        error('You have input an invalid type')
end
        
switch protocol
    case 1 %Parallel type dual superchannels
        protocol='parallel'
        switch k
            case 2
                OMEGA = PermuteSystems(OMEGA,[1 2 4 3 5 6],[d d d d d d]);
            case 3
                OMEGA = PermuteSystems(OMEGA,[1 2 4 6 3 5 7 8],[d d d d d d d d]);
        end
        cvx_begin SDP
            variable W_PIO(d^(2*(k+1)-1),d^(2*(k+1)-1)) semidefinite
             variable rho(d,d) semidefinite
             nofutureDIMS=d*ones(1,2*(k+1)-1);
             PartialTrace(W_PIO,(2+k):2*(k+1)-1,nofutureDIMS) == kron(rho,eye(d^k));
             OMEGA<=kron(W_PIO,eye(d));
             minimise trace(W_PIO)/d^k
        cvx_end
    case 2 %Sequential type dual superchannels
        protocol='sequential'
        switch k
            case 1
             cvx_begin SDP
             variable W_PIO(d^(2*(k+1)-1),d^(2*(k+1)-1)) semidefinite
             variable rho(d,d) semidefinite
             nofutureDIMS=d*ones(1,2*(k+1)-1);
             PartialTrace(W_PIO,(2+k):2*(k+1)-1,nofutureDIMS) == kron(rho,eye(d^k));
             OMEGA<=kron(W_PIO,eye(d));
             minimise trace(W_PIO)/d^k
             cvx_end
            case 2
             cvx_begin SDP
             variable W_PIO(d^(2*(k+1)-1),d^(2*(k+1)-1)) semidefinite
             variable W1(d^3,d^3) semidefinite
             variable rho(d,d) semidefinite
             PartialTrace(W_PIO,5,[d d d d d]) == kron(W1,eye(d));
             PartialTrace(W1,3,[d d d]) == kron(rho,eye(d));
             OMEGA<=kron(W_PIO,eye(d));
             minimise trace(W_PIO)/d^k
             cvx_end
            case 3
             cvx_begin SDP
             variable W_PIO(d^(2*(k+1)-1),d^(2*(k+1)-1)) semidefinite
             variable W1(d^3,d^3) semidefinite
             variable W2(d^5,d^5) semidefinite
             variable rho(d,d) semidefinite
             PartialTrace(W_PIO,7,[d d d d d d d]) == kron(W2,eye(d));
             PartialTrace(W2,5,[d d d d d]) == kron(W1,eye(d));
             PartialTrace(W1,3,[d d d]) == kron(rho,eye(d));
             OMEGA<=kron(W_PIO,eye(d));
             minimise trace(W_PIO)/d^k
             cvx_end
        end
    case 3 %General type superchannels
        protocol='general'
       switch k
            case 1
             cvx_begin SDP
             variable W_PIO(d^(2*(k+1)-1),d^(2*(k+1)-1)) semidefinite
             variable rho(d,d) semidefinite
             nofutureDIMS=d*ones(1,2*(k+1)-1);
             PartialTrace(W_PIO,(2+k):2*(k+1)-1,nofutureDIMS) == kron(rho,eye(d^k));
             OMEGA<=kron(W_PIO,eye(d));
             minimise trace(W_PIO)/d^k
             cvx_end
            case 2
             cvx_begin SDP
             variable W_PIO(d^(2*(k+1)-1),d^(2*(k+1)-1)) semidefinite
             TR(W_PIO,5,[d d d d d]) == TR(W_PIO,[4 5],[d d d d d]);
             TR(W_PIO,3,[d d d d d]) == TR(W_PIO,[2 3],[d d d d d]);
             OMEGA<=kron(W_PIO,eye(d));
             minimise trace(W_PIO)/d^k
             cvx_end
            case 3
             cvx_begin SDP
             variable W_PIO(d^(2*(k+1)-1),d^(2*(k+1)-1)) semidefinite
             TR(W_PIO,7,[d d d d d d d]) == TR(W_PIO,[6 7],[d d d d d d d]);
             TR(W_PIO,5,[d d d d d d d]) == TR(W_PIO,[4 5],[d d d d d d d]);
             TR(W_PIO,3,[d d d d d d d]) == TR(W_PIO,[2 3],[d d d d d d d]);
             OMEGA<=kron(W_PIO,eye(d));
             minimise trace(W_PIO)/d^k
             cvx_end
        end
    otherwise
        error('The type selected does not correspond to any valid superchannel')
end

fprintf('\n\n    The optimal %s unitary %s with d=%i, k=%i, and DUAL method is \n\n',protocol,type,d,k)
F=trace(W_PIO)/d^k
C=W_PIO/trace(W_PIO)*d^k; 

total_time_in_minutes=toc/60
  