classdef G09SCF < handle
    
    properties (Access = private)
        
        verbose;
        info;
        workPath;
        workGjf;
        workLog;
        workDat;
        
        numElecAB;
        overlap;
        toOr;
        harrisMO;
        
        maxSCFIter;
        thresRmsDiffDens;
        thresMaxDiffDens;
        thresDiffEnergy
        
    end
    
    methods (Access = public)
        
        function self = G09SCF(xyz, info, verbose)
            if nargin < 3
                verbose = false;
            end
            self.verbose = verbose;
            self.info = info;
            self.info.xyz = xyz;
            
            self.workPath = tempname();
            system(['mkdir -p ', self.workPath]);
            self.workGjf = [self.workPath, '/run.gjf'];
            self.workLog = [self.workPath, '/run.log'];
            self.workDat = [self.workPath, '/run.dat'];
            self.RunG09('iop(4/199=1) guess=harris');
            
            fLog = fopen(self.workLog, 'r');
            line = fgetl(fLog);
            while ischar(line)
                line = fgetl(fLog);
                if ~isempty(strfind(line, 'alpha electrons'))
                    self.numElecAB = str2double(regexp(line, '\d+', 'match'));
                end
            end
            fclose(fLog);
            
            fDat = fopen(self.workDat, 'rb');
            nbf = sqrt(fread(fDat, 1, 'int') / 8);
            self.overlap = reshape(fread(fDat, nbf^2, 'double'), nbf, nbf);
            fread(fDat, 1, 'int');
            fread(fDat, 1, 'int');
            self.harrisMO = reshape(fread(fDat, nbf^2, 'double'), nbf, nbf);
            fread(fDat, 1, 'int');
            fclose(fDat);
            
            self.toOr = self.ToOrtho(self.overlap);
            
            self.maxSCFIter = 200;
            self.thresRmsDiffDens = 1e-8;
            self.thresMaxDiffDens = 1e-6;
            self.thresDiffEnergy = 1e-6;
        end
        
        function [energy, occOrbList, fockList] = RunSCF(self, guess)
            if nargin < 2
                guess = 'harris';
            end
            if ischar(guess)
                if strcmpi(guess, 'harris')
                    occOrbList = self.GuessHarris();
                end
            elseif iscell(guess)
                occOrbList = guess;
            end
            densList = self.OccOrbToDens(occOrbList);
            lciis = SCF.LCIIS(self.overlap, self.toOr, length(densList), self.verbose);
            energy = 0.0;
            for numIter = 1:self.maxSCFIter
                if self.verbose
                    fprintf('scf iter %f; energy: %f\n', numIter, energy);
                end
                oldEnergy = energy;
                [fockList, energy] = self.FockEnergy(densList);
                fockList = lciis.NewFock(fockList, densList);
                oldDensList = densList;
                occOrbList = self.FockToOccOrb(fockList);
                densList = self.OccOrbToDens(occOrbList);
                if self.Converged(densList, oldDensList, energy, oldEnergy)
                    break;
                end
            end
            if self.verbose
                fprintf('scf done; energy: %f\n', energy);
            end
        end
        
        function delete(self)
            system(['rm -r ', self.workPath]);
        end
        
        function [orbEigVal, orb] = SolveFock(self, fock)
            orFock = self.toOr' * fock * self.toOr;
            [orOrb, orbEigVal] = eig(orFock);
            [orbEigVal, argsort] = sort(diag(orbEigVal));
            orb = self.toOr * orOrb(:, argsort);
        end
        
    end
    
    methods (Access = private)
        
        function RunG09(self, keyword)
            info_ = self.info;
            if isfield(info_, 'numcores')
                numCPUCore = info_.numcores;
            else
                numCPUCore = 1;
            end
            if isfield(info_, 'memory')
                memory = info_.memory;
            else
                memory = '1gb';
            end
            if isfield(info_, 'dft')
                method = info_.dft;
            else
                method = 'dft';
            end
            gjf = {sprintf('%%nprocshared=%d', numCPUCore)};
            gjf = [gjf, {sprintf('%%mem=%s', memory)}];
            gjf = [gjf, {['#p ', method, ' ', info_.basis, ' ', ...
                          'iop(5/13=1, 5/18=-2) scf(maxcycle=1, vshift=-1) ', ...
                          keyword]}];
            gjf = [gjf, {'', 'dummy title', ''}];
            gjf = [gjf, {sprintf('%3d %3d', info_.charge, info_.mult)}];
            for i = 1:size(info_.xyz, 1)
                line = info_.xyz(i, :);
                gjf = [gjf, {sprintf('%3d %15.10f %15.10f %15.10f', line(:))}]; %#ok
            end
            gjf = [gjf, {'', '', ''}];
            gjf = strjoin(gjf, '\n');
            
            fGjf = fopen(self.workGjf, 'w');
            fprintf(fGjf, '%s', gjf);
            fclose(fGjf);
            
            system(['g09binfile=', self.workDat, ' ' ...
                    'g09 ', self.workGjf, ' ', self.workLog]);
            
            fLog = fopen(self.workLog, 'r');
            line = fgetl(fLog);
            while ischar(line)
                prevLine = line;
                line = fgetl(fLog);
            end
            fclose(fLog);
            if isempty(strfind(prevLine, 'Normal termination'))
                throw(MException('G09SCF:g09failed', 'g09 failed'));
            end
        end
        
        function toOr = ToOrtho(~, overlap)
            [eigVec, eigVal] = eig(overlap);
            eigVal = reshape(diag(eigVal), 1, []);
            keep = eigVal > 1.0e-6;
            toOr = eigVec(:, keep) ./ repmat(sqrt(eigVal(keep)), size(eigVec, 1), 1);
        end
        
        function occOrbList = GuessHarris(self)
            occOrbList = cell(1, length(unique(self.numElecAB)));
            for i = 1:length(unique(self.numElecAB))
                occOrbList{i} = self.harrisMO(:, 1:self.numElecAB(i));
            end
        end
        
        function [fockList, energy] = FockEnergy(self, densList)
            nbf = size(self.overlap, 1);
            fDat = fopen(self.workDat, 'wb');
            for densCell = densList
                dens = densCell{1};
                fwrite(fDat, 8 * nbf^2, 'int');
                fwrite(fDat, dens, 'double');
                fwrite(fDat, 8 * nbf^2, 'int');
            end
            fclose(fDat);
            self.RunG09('iop(5/199=1) guess=core')
            fDat = fopen(self.workDat, 'rb');
            fread(fDat, 1, 'int');
            energy = fread(fDat, 1, 'double');
            fread(fDat, 1, 'int');
            fockList = cell(1, length(densList));
            for i = 1:length(densList)
                fread(fDat, 1, 'int');
                fockList{i} = reshape(fread(fDat, nbf^2, 'double'), nbf, nbf);
                fread(fDat, 1, 'int');
            end
            fclose(fDat);
        end
        
        function occOrbList = FockToOccOrb(self, fockList)
            occOrbList = cell(1, length(fockList));
            for i = 1:length(fockList)
                [~, orb] = self.SolveFock(fockList{i});
                occOrbList{i} = orb(:, 1:self.numElecAB(i));
            end
        end
        
        function densList = OccOrbToDens(~, occOrbList)
            densList = cell(1, length(occOrbList));
            for i = 1:length(occOrbList)
                densList{i} = occOrbList{i} * occOrbList{i}';
            end
        end
        
        function con = Converged(self, densList, oldDensList, energy, oldEnergy)
            diffDens = zeros(numel(densList{1}), length(densList));
            for i = 1:length(densList)
                diffDens(:, i) = reshape(densList{i} - oldDensList{i}, [], 1);
            end
            diffEnergy = abs(energy - oldEnergy);
            rmsDiffDens = sqrt(mean(mean(diffDens.^2)));
            maxDiffDens = max(max(abs(diffDens)));
            if self.verbose
                form = '  %s: %.3e thres: %.3e\n';
                fprintf(form, 'rmsDiffDens', rmsDiffDens, self.thresRmsDiffDens);
                fprintf(form, 'maxDiffDens', maxDiffDens, self.thresMaxDiffDens);
                fprintf(form, 'diffEnergy ', diffEnergy, self.thresDiffEnergy);
            end
            con = rmsDiffDens < self.thresRmsDiffDens && ...
                  maxDiffDens < self.thresMaxDiffDens && ...
                   diffEnergy < self.thresDiffEnergy;
        end
        
    end
    
end