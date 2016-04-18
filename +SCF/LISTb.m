classdef LISTb < handle
    
    properties (Access = private)
        
        fockInList;
        
        packs;
        maxNumPack;
        
        verbose;
        
    end
    
    methods
        
        function self = LISTb(verbose, maxNumPack)
            if(nargin < 2)
                maxNumPack = 5;
            end
            self.packs = {};
            self.maxNumPack = maxNumPack;
            self.verbose = verbose;
        end
        
        function newFockList = NewFock(self, fockList, densList, energy)
            newPack.fockList = fockList;
            newPack.energy = energy;
            newPack.densList = densList;
            
            if ~isempty(self.packs)
                dFVec = [];
                for i = 1:length(fockList)
                    temp = reshape(fockList{i} - self.fockInList{i}, [], 1);
                    dFVec = [dFVec; temp]; %#ok
                end
                newPack.dFVec = dFVec / 2;
            end
            
            if length(self.packs) == self.maxNumPack
                self.packs(1:(end - 1)) = self.packs(2:end);
                self.packs{end} = newPack;
            else
                self.packs{end + 1} = newPack;
            end
            
            usePack = self.packs(2:end);
            numPack = length(usePack);
            if numPack <= 1
                newFockList = fockList;
                self.fockInList = newFockList;
                return;
            end
            
            hess = zeros(numPack, numPack);
            for i = 1:numPack
                for j = 1:numPack
                    densi = usePack{i}.densList;
                    densj = usePack{j}.densList;
                    dDij = [];
                    for s = 1:length(fockList)
                        dDij = [dDij; reshape(densi{s} - densj{s}, [], 1)]; %#ok
                    end
                    hess(i, j) = usePack{i}.energy + dDij' * usePack{i}.dFVec;
                end
            end
            onesVec = ones(numPack, 1);
            hessL = [hess', onesVec; onesVec', 0];
            [coeff, ~] = linsolve(hessL, [zeros(numPack, 1); 1]);
            coeff = coeff(1:(end - 1));
            
            if self.verbose
                fprintf('  coeff:')
                fprintf(' %6.3f', coeff);
                fprintf('\n')
            end
            
            newFockList = cell(1, length(fockList));
            for spin = 1:length(fockList)
                allFock = zeros(numel(fockList{1}), numPack);
                for useInd = 1:numPack
                    fock = usePack{useInd}.fockList{spin};
                    allFock(:, useInd) = reshape(fock, [], 1);
                end
                newFockList{spin} = reshape(allFock * coeff, size(fockList{1}));
            end
            
            self.fockInList = newFockList;
            
        end
        
    end
    
end
