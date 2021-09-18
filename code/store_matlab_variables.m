%Script which creates stores all variables used for the computer assisted proof
clear
also_consider_qutrits=0; %Be careful, qutrit calculations consume a big ammount of Ram memory (more than 10 GB) and a long time (maybe more than 1h)

%For qubits, we only need to consider unitary transposition (complex
%conjugation can be done for free)

% k=1 for qubits
d=2;k=1;type=1;protocol=1;dual=0;
[F,Sd2k1type1protocol1] = optimal_fU(d,k,protocol,type,dual);
dual=1;
[F,Wd2k1type1protocol1] = optimal_fU(d,k,protocol,type,dual);

% k=2 for qubits
d=2;k=2;type=1;protocol=1;dual=0;
[F,Sd2k2type1protocol1] = optimal_fU(d,k,protocol,type,dual);
d=2;k=2;type=1;protocol=1;dual=1;
[F,Wd2k2type1protocol1] = optimal_fU(d,k,protocol,type,dual);

d=2;k=2;type=1;protocol=2;dual=0;
[F,Sd2k2type1protocol2] = optimal_fU(d,k,protocol,type,dual);
dual=1;
[F,Wd2k2type1protocol2] = optimal_fU(d,k,protocol,type,dual);

d=2;k=2;type=1;protocol=3;dual=0;
[F,Sd2k2type1protocol3] = optimal_fU(d,k,protocol,type,dual);
dual=1;
[F,Wd2k2type1protocol3] = optimal_fU(d,k,protocol,type,dual);

mkdir MatlabVariables %Creates the folder there the variables will be saved
cd MatlabVariables/   %Accesses the folder there the variables will be saved

 save 'Sd2k1type1protocol1' Sd2k1type1protocol1
 save 'Wd2k1type1protocol1' Wd2k1type1protocol1
 
 save 'Sd2k2type1protocol1' Sd2k2type1protocol1
 save 'Wd2k2type1protocol1' Wd2k2type1protocol1
 
 save 'Sd2k2type1protocol2' Sd2k2type1protocol2
 save 'Wd2k2type1protocol2' Wd2k2type1protocol2
 
 save 'Sd2k2type1protocol3' Sd2k2type1protocol3
 save 'Wd2k2type1protocol3' Wd2k2type1protocol3
 cd .. % Returns to the original folder
 
 disp('The qubit variables were saved in the folder MatlabVariables');

if also_consider_qutrits
d=3;k=1;type=1;protocol=1;dual=0;
[F,Sd3k1type1protocol1] = optimal_fU(d,k,protocol,type,dual);
dual=1;
[F,Wd3k1type1protocol1] = optimal_fU(d,k,protocol,type,dual);

d=3;k=2;type=1;protocol=1;dual=0;
[F,Sd3k2type1protocol1] = optimal_fU(d,k,protocol,type,dual);
dual=1;
[F,Wd3k2type1protocol1] = optimal_fU(d,k,protocol,type,dual);

d=3;k=2;type=1;protocol=2;dual=0;
[F,Sd3k2type1protocol2] = optimal_fU(d,k,protocol,type,dual);
dual=1;
[F,Wd3k2type1protocol2] = optimal_fU(d,k,protocol,type,dual);

d=3;k=2;type=1;protocol=3;dual=0;
[F,Sd3k2type1protocol3] = optimal_fU(d,k,protocol,type,dual);
dual=1;
[F,Wd3k2type1protocol3] = optimal_fU(d,k,protocol,type,dual);

 cd MatlabVariables/   %Accesses the folder there the variables will be saved

 save 'Sd3k1type1protocol1' Sd3k1type1protocol1
 save 'Wd3k1type1protocol1' Wd3k1type1protocol1
 
 save 'Sd3k2type1protocol1' Sd3k2type1protocol1
 save 'Wd3k2type1protocol1' Wd3k2type1protocol1
 
 save 'Sd3k2type1protocol2' Sd3k2type1protocol2
 save 'Wd3k2type1protocol2' Wd3k2type1protocol2
 
 save 'Sd3k2type1protocol3' Sd3k2type1protocol3
 save 'Wd3k2type1protocol3' Wd3k2type1protocol3
 
 cd .. % Returns to the original folder
 disp('The qutrit unitary transposition variables were saved in the folder MatlabVariables');

% For unitary Inversion

d=3;k=1;type=3;protocol=1;dual=0;
[F,Sd3k1type3protocol1] = optimal_fU(d,k,protocol,type,dual);
dual=1;
[F,Wd3k1type3protocol1] = optimal_fU(d,k,protocol,type,dual);
 

d=3;k=2;type=3;protocol=2;dual=0;
[F,Sd3k2type3protocol2] = optimal_fU(d,k,protocol,type,dual);
dual=1;
[F,Wd3k2type3protocol2] = optimal_fU(d,k,protocol,type,dual);

d=3;k=2;type=3;protocol=3;dual=0;
[F,Sd3k2type3protocol3] = optimal_fU(d,k,protocol,type,dual);
dual=1;
[F,Wd3k2type3protocol3] = optimal_fU(d,k,protocol,type,dual);


 cd MatlabVariables/   %Accesses the folder there the variables will be saved
 
 save 'Sd3k1type3protocol1' Sd3k1type3protocol1
 save 'Wd3k1type3protocol1' Wd3k1type3protocol1

 save 'Sd3k2type3protocol1' Sd3k2type3protocol1
 save 'Wd3k2type3protocol1' Wd3k2type3protocol1
 
 save 'Sd3k2type3protocol2' Sd3k2type3protocol2
 save 'Wd3k2type3protocol2' Wd3k2type3protocol2
 
 save 'Sd3k2type3protocol3' Sd3k2type3protocol3
 save 'Wd3k2type3protocol3' Wd3k2type3protocol3
 cd .. % Returns to the original folder
 disp('The qutrit unitary inversion variables were saved in the folder MatlabVariables');
end