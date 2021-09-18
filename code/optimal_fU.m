function [F,C,OMEGA] = optimal_fU(d,k,protocol,type,dual)
%Function that evaluates the maximal average fidelity for unitary transformations.
%protocol: 1 for parallel, 2 for sequential, 3 for general
%typoe: %1 for Transposition, 2 for Conjugate, 3 for Inversion 
%dual: 0 for primal and 1 for dual

OMEGA= make_OMEGA_permutation(d,k,type);
switch type
    case 1
        type='transposition';
    case 2
        type='conjugate';
    case 3
        type='inversion';
    otherwise
        error('You have input an invalid type')
end

if dual==0  %Start the code for primal
    reorderSuperchannel=0;
switch protocol
    case 1 %Parallel type superchannels
        reorderSuperchannel=1;
        protocol='parallel';
        cvx_begin SDP quiet
            variable C(d^(2*(k+1)),d^(2*(k+1))) semidefinite
             switch k
             case 2
                OMEGA = PermuteSystems(OMEGA,[1 2 4 3 5 6],[d d d d d d]);
             case 3
                OMEGA = PermuteSystems(OMEGA,[1 2 4 6 3 5 7 8],[d d d d d d d d]);
             end
             PartialTrace(C,4,[d d^k d^k d]) == kron(PartialTrace(C,[3 4],[d d^k d^k d]),eye(d^k)/d^k);
             PartialTrace(C,[2 3 4],[d d^k d^k d])  == eye(d)*d^k;
             maximise real(OMEGA(:)'*C(:))
        cvx_end    
    case 2 %Sequential type superchannels
        cvx_begin SDP quiet quiet
            protocol='sequential';
            variable C(d^(2*(k+1)),d^(2*(k+1))) semidefinite
            switch k
                case 1
                     PartialTrace(C,4,[d d d d])==kron(PartialTrace(C,[3 4],[d d d d]),eye(d)/d);
                     PartialTrace(C,[2 3 4],[d d d d])==eye(d)*d;
                case 2
                     PartialTrace(C,6,[d d d d d d]) == kron(PartialTrace(C,[5 6],[d d d d d d]),eye(d)/d);
                     PartialTrace(C,[4 5 6],[d d d d d d]) == kron(PartialTrace(C,[3 4 5 6],[d d d d d d]),eye(d)/d);
                     PartialTrace(C,[2 3 4 5 6],[d d d d d d])==eye(d)*d^2;
                case 3
                     PartialTrace(C,8, [d d d d d d d d])==kron(PartialTrace(C,[7 8], [d d d d d d d d]),eye(d)/d);
                     PartialTrace(C,[6 7 8], [d d d d d d d d])==kron(PartialTrace(C,[5 6 7 8], [d d d d d d d d]),eye(d)/d);
                     PartialTrace(C,[4 5 6 7 8], [d d d d d d d d])==kron(PartialTrace(C,[3 4 5 6 7 8], [d d d d d d d d]),eye(d)/d);
                     PartialTrace(C,[2 3 4 5 6 7 8], [d d d d d d d d])==eye(d)*d^3;     
            end
            maximise real(OMEGA(:)'*C(:))
        cvx_end
    case 3 %General type superchannels
        cvx_begin SDP quiet
            protocol='general';
            variable C(d^(2*(k+1)),d^(2*(k+1))) semidefinite
            is_general_superchannel(C,d,k)
            maximise real(OMEGA(:)'*C(:))
        cvx_end
    otherwise
        error('The type selected does not correspond to any valid superchannel')
end

fprintf('\n\n    The optimal %s unitary %s with d=%i and k=%i is \n\n',protocol,type,d,k)
F=OMEGA(:)'*C(:)

if reorderSuperchannel
             switch k
             case 2
                C = PermuteSystems(C,[1 2 4 3 5 6],[d d d d d d],0,1);
             case 3
                C = PermuteSystems(C,[1 2 4 6 3 5 7 8],[d d d d d d d d],0,1);
             end
end

else  %write the dual part
    reorderSuperchannel=0;
    switch protocol
    case 1 %Parallel type dual superchannels
        protocol='parallel';
        reorderSuperchannel=1;
        switch k
            case 2
                OMEGA = PermuteSystems(OMEGA,[1 2 4 3 5 6],[d d d d d d]);
            case 3
                OMEGA = PermuteSystems(OMEGA,[1 2 4 6 3 5 7 8],[d d d d d d d d]);
        end
        cvx_begin SDP quiet
            variable W_PIO(d^(2*(k+1)-1),d^(2*(k+1)-1)) semidefinite
             variable rho(d,d) semidefinite
             nofutureDIMS=d*ones(1,2*(k+1)-1);
             PartialTrace(W_PIO,(2+k):2*(k+1)-1,nofutureDIMS) == kron(rho,eye(d^k));
             OMEGA<=kron(W_PIO,eye(d));
             minimise trace(W_PIO)/d^k
        cvx_end
    case 2 %Sequential type dual superchannels
        protocol='sequential';
        switch k
            case 1
             cvx_begin SDP quiet
             variable W_PIO(d^(2*(k+1)-1),d^(2*(k+1)-1)) semidefinite
             variable rho(d,d) semidefinite
             nofutureDIMS=d*ones(1,2*(k+1)-1);
             PartialTrace(W_PIO,(2+k):2*(k+1)-1,nofutureDIMS) == kron(rho,eye(d^k));
             OMEGA<=kron(W_PIO,eye(d));
             minimise trace(W_PIO)/d^k
             cvx_end
            case 2
             cvx_begin SDP quiet
             variable W_PIO(d^(2*(k+1)-1),d^(2*(k+1)-1)) semidefinite
             variable W1(d^3,d^3) semidefinite
             variable rho(d,d) semidefinite
             PartialTrace(W_PIO,5,[d d d d d]) == kron(W1,eye(d));
             PartialTrace(W1,3,[d d d]) == kron(rho,eye(d));
             OMEGA<=kron(W_PIO,eye(d));
             minimise trace(W_PIO)/d^k
             cvx_end
            case 3
             cvx_begin SDP quiet
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
        protocol='general';
       switch k
            case 1
             cvx_begin SDP quiet
             variable W_PIO(d^(2*(k+1)-1),d^(2*(k+1)-1)) semidefinite
             variable rho(d,d) semidefinite
             nofutureDIMS=d*ones(1,2*(k+1)-1);
             PartialTrace(W_PIO,(2+k):2*(k+1)-1,nofutureDIMS) == kron(rho,eye(d^k));
             OMEGA<=kron(W_PIO,eye(d));
             minimise trace(W_PIO)/d^k
             cvx_end
            case 2
             cvx_begin SDP quiet
             variable W_PIO(d^(2*(k+1)-1),d^(2*(k+1)-1)) semidefinite
             TR(W_PIO,5,[d d d d d]) == TR(W_PIO,[4 5],[d d d d d]);
             TR(W_PIO,3,[d d d d d]) == TR(W_PIO,[2 3],[d d d d d]);
             OMEGA<=kron(W_PIO,eye(d));
             minimise trace(W_PIO)/d^k
             cvx_end
            case 3
             cvx_begin SDP quiet
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
C=W_PIO/trace(W_PIO)*d^k; %Assing the dual affine superchannel to the variable C

if reorderSuperchannel
             switch k
             case 2
                C = PermuteSystems(C,[1 2 4 3 5],[d d d d d],0,1);
             case 3
                C = PermuteSystems(C,[1 2 4 6 3 5 7],[d d d d d d d],0,1);
             end
end


end %end else
end %end function