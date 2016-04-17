classdef ADIIS < SCF.EDIIS
    
    methods (Access = public)
        
        function self = ADIIS(verbose, maxNumPack)
            self@SCF.EDIIS(verbose, maxNumPack);
        end
        
    end
    
    methods (Access = protected)
        
        function [hessian, linPart] = HessianAndLinPart(self)
            numPack = length(self.packs);
            numSpin = length(self.packs{1}.fockList);
            errFVec = cell(1, numSpin);
            errDVec = cell(1, numSpin);
            nbfsq = numel(self.packs{1}.fockList{1});
            lastPack = self.packs{end};
            
            % compute errors wrt. the latest Fock or density
            for s = 1:numSpin
                errFVec{s} = zeros(nbfsq, numPack);
                errDVec{s} = zeros(nbfsq, numPack);
                for i = 1:numPack
                    errF = self.packs{i}.fockList{s} - lastPack.fockList{s};
                    errD = self.packs{i}.densList{s} - lastPack.densList{s};
                    errFVec{s}(:, i) = reshape(errF, [], 1);
                    errDVec{s}(:, i) = reshape(errD, [], 1);
                end
            end
            
            linPart = zeros(numPack, 1);
            hessian = zeros(numPack, numPack);
            for s = 1:numSpin
                lastFockVec = reshape(lastPack.fockList{s}, [], 1);
                linPart = linPart + 2 * (lastFockVec' * errDVec{s})';
                hessian = hessian + errFVec{s}' * errDVec{s};
            end
            hessian = hessian + hessian';
        end
        
    end
    
end