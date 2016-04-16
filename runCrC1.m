cart = [...
    24                0.00000000    0.00000000    0.40000000
    6                 0.00000000    0.00000000   -1.60000000];

import SCF.*;


info.charge = 0;
info.multiplicity = 3;
info.cartesian = cart;
info.basisSet = '6-31g';
info.dft = 'b3lyp';
diisType = 'L20';

scf = UKS(info);

% scf.RunG09();
[guessDensVec, guessOrbital] = scf.CoreGuess();
guess.guessOrbital = guessOrbital;

% guessDensVec = scf.SADGuess();
% guess.guessDensVec = [guessDensVec];

output = scf.SCF(guess, diisType);
output

