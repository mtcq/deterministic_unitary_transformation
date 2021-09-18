clear 
tic
%This is a script that evaluates the maximal average fidelity for unitary transformations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Start of Adjustable Parameters  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d=2               %Dimension 'd'
        k=2             %Number of uses 'k'
        protocol=2; %1 for parallel, 2 for sequential, 3 for general
        type=1;      %1 for Transposition, 2 for Conjugate, 3 for Inversion 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   End of Adjustable Parameters   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
OMEGA= make_OMEGA_permutation(d,k,type);
%OMEGA=make_OMEGA_group(d,k,type);


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

reorderSuperchannel=0; %Variable that indicates if we should reoder the superchannel C (useful for parallel superchannels)
switch protocol
    case 1 %Parallel type superchannels
        reorderSuperchannel=1;
        protocol='parallel'
        cvx_begin SDP
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
        cvx_begin SDP
            protocol='sequential'
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
        cvx_begin SDP
            protocol='general'
            variable C(d^(2*(k+1)),d^(2*(k+1))) semidefinite
            is_general_superchannel(C,d,k)
            maximise real(OMEGA(:)'*C(:))
        cvx_end
    otherwise
        error('The type selected does not correspond to any valid superchannel')
end

fprintf('\n\n    The optimal %s unitary %s with d=%i and k=%i is \n\n',protocol,type,d,k)
%F=trace(OMEGA*C)
F=OMEGA(:)'*C(:)
if reorderSuperchannel
             switch k
             case 2
                C = PermuteSystems(C,[1 2 4 3 5 6],[d d d d d d],0,1);
             case 3
                C = PermuteSystems(C,[1 2 4 6 3 5 7 8],[d d d d d d d d],0,1);
             end
end


total_time_in_minutes=toc/60
  