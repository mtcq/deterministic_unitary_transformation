function OMEGA=make_OMEGA_permutation(d,k,type)
%This function makes the performance operator omega using the basis
%constructe via group representation methods. It only works for k<=3, but
%it can easily be adapted for any k
%type=1 is TRANS, 2 is CONJ, 3 is INV

% Omega will be outputed in the order: P I_1 O_1 I_2 O_2 F
%Paper notation:   P   I_1   O_1   I_2   O_2   F
%Matlab notation:  1    2     3     4     5    6

if k>=4
    error('The permutation method for generating OMEGA only holds for k=1, k=2, and k=3')
    pause
end

if type==1
    %'type=1, Make Omega for Unitary Transposition!'
    OMEGA=0;
 switch k
     case 1
        P=basis_Ukstar_overcomplete_method(d,k+1,[1]);
        for i=1:size(P,3)
          OMEGA = OMEGA +  kron(P(:,:,i),conj(P(:,:,i)));
        end
          OMEGA = PermuteSystems(OMEGA,[3 1 2 4],[d d d d],0,1)/d^2;
          
     case 2
          P=basis_Ukstar_overcomplete_method(d,k+1,[1 2]);
        for i=1:size(P,3)
          OMEGA = OMEGA +  kron(P(:,:,i),conj(P(:,:,i)));
        end
          OMEGA = PermuteSystems(OMEGA,[3 5 1 2 4 6],[d d d d d d],0,1)/d^2;
     case 3
          P=basis_Ukstar_overcomplete_method(d,k+1,[1 2 3]);
        for i=1:size(P,3)
          OMEGA = OMEGA +  kron(P(:,:,i),conj(P(:,:,i)));
        end
           OMEGA = PermuteSystems(OMEGA,[3 5 7 1 2 4 6 8],[d d d d d d d d],0,1)/d^2;
 end
 
end

if type==2
    %'type=2, Make Omega for Unitary Conjugation!'

      OMEGA=0;
        P=basis_Uk_overcomplete_method(d,k+1);
        for i=1:size(P,3)
            if trace(P(:,:,i)*P(:,:,i)')>0.0001                     
                OMEGA = OMEGA +  kron(conj(P(:,:,i)),P(:,:,i))/trace(P(:,:,i)*P(:,:,i)');
            end
        end
 switch k
     case 1
          OMEGA = PermuteSystems(OMEGA,[2 1 3 4],[d d d d],0,1)/d^2;
     case 2
          OMEGA = PermuteSystems(OMEGA,[2 4 1 3 5 6],[d d d d d d],0,1)/d^2;  %Use 1 to invert 
     case 3
          OMEGA = PermuteSystems(OMEGA,[2 4 6 1 3 5 7 8],[d d d d d d d d],0,1)/d^2;
 end
end

if type==3
    %'type=3, Make Omega for Unitary Inversion!'
        OMEGA=0;
        P=basis_Uk_overcomplete_method(d,k+1);
        for i=1:size(P,3)
            if trace(P(:,:,i)'*P(:,:,i))>0.0001                     
                OMEGA = OMEGA + kron(conj(P(:,:,i)),P(:,:,i))/trace(P(:,:,i)*P(:,:,i)');
            end
        end
        
 switch k
     case 1
         OMEGA = PermuteSystems(OMEGA,[3 1 2 4],[d d d d],0,1)/d^2;
     case 2
         %OMEGA = PermuteSystems(OMEGA,[1 5 2 6 3 4],[d d d d d d])/d^2;
         OMEGA = PermuteSystems(OMEGA,[3 5 1 2 4 6],[d d d d d d],0,1)/d^2;
     case 3
         OMEGA = PermuteSystems(OMEGA,[3 5 7 1 2 4 6 8],[d d d d d d d d],0,1)/d^2;
 end
end
end









