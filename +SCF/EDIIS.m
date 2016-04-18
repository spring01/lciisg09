classdef EDIIS < handle
    
    properties (Access = protected)
        
        packs;
        
    end
    
    properties (Access = private)
        
        verbose;
        maxNumPack;
        
    end
    
    methods (Access = public)
        
        function self = EDIIS(verbose, maxNumPack)
            if(nargin < 1)
                verbose = false;
            end
            if(nargin < 2)
                maxNumPack = 20;
            end
            self.packs = {};
            self.verbose = verbose;
            self.maxNumPack = maxNumPack;
        end
        
        function newFockList = NewFock(self, fockList, densList, energy)
            newPack.fockList = fockList;
            newPack.densList = densList;
            newPack.energy = energy;
            if length(self.packs) == self.maxNumPack
                self.packs(1:(end - 1)) = self.packs(2:end);
                self.packs{end} = newPack;
            else
                self.packs{end + 1} = newPack;
            end
            
            numPack = length(self.packs);
            if numPack == 1
                newFockList = fockList;
                return;
            end
            
            [hessian, linPart] = self.HessianAndLinPart();
            
            options = optimoptions('quadprog', 'Algorithm', 'active-set', ...
                                   'Display', 'off');
            coeff = quadprog(hessian, linPart, -eye(numPack), ...
                             zeros(numPack, 1),  ones(1, numPack), 1, ...
                             [], [], [], options);
            if self.verbose
                fprintf('  coeff:');
                fprintf(' %6.3f', coeff);
                fprintf('\n')
            end
            
            newFockList = cell(1, length(fockList));
            for spin = 1:length(fockList)
                allFock = zeros(numel(fockList{1}), numPack);
                for useInd = 1:numPack
                    fock = self.packs{useInd}.fockList{spin};
                    allFock(:, useInd) = reshape(fock, [], 1);
                end
                newFockList{spin} = reshape(allFock * coeff, size(fockList{1}));
            end
        end
        
    end
    
    methods (Access = protected)
        
        function [hessian, linPart] = HessianAndLinPart(self)
            numPack = length(self.packs);
            numSpin = length(self.packs{1}.fockList);
            hessian = zeros(numPack, numPack);
            linPart = zeros(numPack, 1);
            for i = 1:numPack
                iFL = self.packs{i}.fockList;
                iDL = self.packs{i}.densList;
                linPart(i) = self.packs{i}.energy;
                for j = 1:numPack
                    jFL = self.packs{j}.fockList;
                    jDL = self.packs{j}.densList;
                    errF = [];
                    errD = [];
                    for sp = 1:length(iFL)
                        errF = [errF; reshape(iFL{sp} - jFL{sp}, [], 1)]; %#ok
                        errD = [errD; reshape(iDL{sp} - jDL{sp}, [], 1)]; %#ok
                    end
                    hessian(j, i) = errF' * errD;
                end
            end
            hessian = -hessian / numSpin;
        end
        
    end
    
end