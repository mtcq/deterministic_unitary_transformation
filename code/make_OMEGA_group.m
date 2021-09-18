function OMEGA=make_OMEGA_group(d,k,type)
%This function makes the performance operator omega using the basis
%constructe via group representation methods. It only works for k<=2
%type=1 is TRANS, 2 is CONJ, 3 is INV

% Omega will be outputed in the order: P I_1 O_1 I_2 O_2 F
%Paper notation:   P   I_1   O_1   I_2   O_2   F
%Matlab notation:  1    2     3     4     5    6

if k>=3
    error('The analytic method for generating OMEGA only holds for k=1 and k=2')
    pause
end

if type==1
    'type=1, Make Omega for Unitary Transposition!'
    if k==1
        mask=[d d d d];
        P1=IsotropicState(d,1);
        P2=eye(d^2)-P1;
        OMEGA = 1/(d^2) * ( kron(P1,conj(P1))/1 + kron(P2,conj(P2))/(d^2-1) );
        %Here OMEGA has the order P O F I (in MATLAB: 1 3 4 2)
        OMEGA = PermuteSystems(OMEGA,[1 3 4 2],mask,0,1);
        %Permute OMEGA using the invert PermuteSystem
    end    
    if k==2
        OMEGA=zeros(d^6);
        %P=basis_UUUstar_PAVIA(d);
        P=basis_UstarUstarU(d);
        for i=1:6
            if trace(P(:,:,i)*P(:,:,i)')>0.0001
                             
                OMEGA = OMEGA +  kron(P(:,:,i),conj(P(:,:,i)))/trace(P(:,:,i)*P(:,:,i)');
            end
        end
     %Here OMEGA has the order P O_1 O_2 F I_1 I_2 (in MATLAB: 1 3 5 6 2 4)
    % OMEGA = PermuteSystems(OMEGA,[1 3 5 6 2 4],[d d d d d d],0,1)/d^2;
     OMEGA = PermuteSystems(OMEGA,[3 5 1 2 4 6],[d d d d d d],0,1)/d^2;
     %Permute OMEGA using the invert PermuteSystem
    end
end

if type==2
    'type=2, Make Omega for Unitary Conjugation!'
    if k==1
        P1=SymmetricProjection(d);
        P2=AntisymmetricProjection(d);
        OMEGA = 1/(d^2) * ( kron(P1,conj(P1))/trace(P1) + kron(P2,conj(P2))/trace(P2) );
        OMEGA = PermuteSystems(OMEGA,[2 1 3 4],[d d d d],0,1);
    end
    if k==2
         OMEGA=zeros(d^6);
        P=basis_UUU(d);
        for i=1:6
            if trace(P(:,:,i)*P(:,:,i)')>0.0001                
                OMEGA = OMEGA +  kron(conj(P(:,:,i)),P(:,:,i))/trace(P(:,:,i)*P(:,:,i)');
            end
        end
        OMEGA = PermuteSystems(OMEGA,[2 4 1 3 5 6],[d d d d d d],0,1)/d^2;  %Use 1 to invert  
    end
end


if type==3
    'type=3, Make Omega for Unitary Inversion!'
    if k==1
        mask=[d d d d];
        P1=SymmetricProjection(d);
        P2=AntisymmetricProjection(d);
        OMEGA = 1/(d^2) * ( kron(P1,conj(P1))/trace(P1) + kron(P2,conj(P2))/trace(P2) );
        OMEGA = PermuteSystems(OMEGA,[3 1 2 4],mask);
    end
    if k==2
        OMEGA=zeros(d^6);
        P=basis_UUU(d);
        for i=1:6
            if trace(P(:,:,i)*P(:,:,i)')>0.0001                     
                OMEGA = OMEGA +  kron(conj(P(:,:,i)),P(:,:,i))/trace(P(:,:,i)*P(:,:,i)');
            end
        end
 % OMEGA = PermuteSystems(OMEGA,[1 5 2 6 3 4],[d d d d d d])/d^2;
  OMEGA = PermuteSystems(OMEGA,[3 5 1 2 4 6],[d d d d d d],0,1)/d^2;
    end
end

end