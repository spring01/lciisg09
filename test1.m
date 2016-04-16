xyz = [8, 0.0,  0.000000,  0.110200
       1, 0.0,  0.711600, -0.440800
       1, 0.0, -0.711600, -0.440800];

info.dft = 'b3lyp';
info.basis = 'sto-3g';
info.charge = 0;
info.mult = 1;

scf = SCF.G09SCF(xyz, info, true);
[energy, occOrbList, fockList] = scf.RunSCF()