classdef LCIIS < handle
    
    properties (Access = private)
        
        packs;
        comms;
        
        toOr;
        orTo;
        
        verbose;
        maxNumPack
        bigMat;
        
        maxIterNewton;
        gradNormThres = 1e-12;
        
    end
    
    methods
        
        function self = LCIIS(overlap, toOr, numSpin, verbose, maxNumPack, ...
                              maxIterNewton)
            if nargin < 4
                verbose = false;
            end
            if nargin < 5
                maxNumPack = 60;
            end
            if nargin < 6
                maxIterNewton = 200;
            end
            self.verbose = verbose;
            self.packs = {};
            self.maxNumPack = maxNumPack;
            self.maxIterNewton = maxIterNewton;
            nbf = size(overlap, 1);
            self.comms = zeros(numSpin * nbf * (nbf - 1) / 2, maxNumPack^2);
            self.bigMat = zeros(maxNumPack^2, maxNumPack^2);
            self.toOr = toOr;
            self.orTo = overlap * toOr;
        end
        
        function newFockList = NewFock(self, fockList, densList, ~)
            if length(self.packs) == self.maxNumPack
                self.PreUpdateFull();
            else
                self.PreUpdateNotFull();
            end
            newPack.fockList = fockList;
            newPack.trFList = cell(1, length(fockList));
            newPack.trDList = cell(1, length(fockList));
            for spin = 1:length(fockList)
                newPack.trFList{spin} = self.toOr' * fockList{spin};
                newPack.trDList{spin} = densList{spin} * self.orTo;
            end
            if length(self.packs) == self.maxNumPack
                self.packs(1:(end - 1)) = self.packs(2:end);
                self.packs{end} = newPack;
            else
                self.packs{end + 1} = newPack;
            end
            self.UpdateComms();
            self.UpdateBigMat();
            
            numPack = length(self.packs);
            shapeT = [numPack, numPack, numPack, numPack];
            tensor = reshape(self.bigMat(1:numPack^2, 1:numPack^2), shapeT);
            tensorGrad = tensor + permute(tensor, [2 1 3 4]);
            tensorHess = tensor + permute(tensor, [1 3 2 4]) ...
                                + permute(tensor, [1 4 2 3]) ...
                                + permute(tensor, [2 1 3 4]) ...
                                + permute(tensor, [2 3 1 4]) ...
                                + permute(tensor, [2 4 1 3]);
            coeff = zeros(numPack, 1);
            for useStart = 1:numPack
                useInd = useStart:numPack;
                tensorUse = tensor(useInd, useInd, useInd, useInd);
                tGrad = tensorGrad(useInd, useInd, useInd, useInd);
                tHess = tensorHess(useInd, useInd, useInd, useInd);
                [ok, iniCoeffUse] = self.InitialCoeffUse(tensorUse);
                if ~ok
                    disp('InitialCoeffUse failed; reducing tensor size');
                    continue;
                end
                [ok, coeffUse] = self.NewtonSolver(tGrad, tHess, iniCoeffUse);
                if ~ok
                    disp('NewtonSolver failed; reducing tensor size');
                    continue;
                else
                    coeff(useInd) = coeffUse(:);
                    break
                end
            end
            if self.verbose
                fprintf('  coeff:')
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
        
        function comm = CommBetween(self, i, j)
            comm = [];
            numSpin = length(self.packs{1}.fockList);
            for spin = 1:numSpin
                fd = self.packs{i}.trFList{spin} * self.packs{j}.trDList{spin};
                commMat = fd - fd';
                comm = [comm; commMat(triu(true(size(commMat)), 1))]; %#ok
            end
        end
        
    end
    
    methods (Access = private)
        
        function PreUpdateFull(self)
            numPack = self.maxNumPack;
            for i = 2:numPack
                srcInd = ((i - 1) * numPack + 2):(i * numPack);
                tgtInd = srcInd - numPack - 1;
                self.comms(:, tgtInd) = self.comms(:, srcInd);
                self.bigMat(tgtInd, :) = self.bigMat(srcInd, :);
                self.bigMat(:, tgtInd) = self.bigMat(:, srcInd);
            end
        end
        
        function PreUpdateNotFull(self)
            numPack = length(self.packs) + 1;
            for i = (numPack - 1):-1:2
                srcInd = ((i - 1) * (numPack - 1) + 1):(i * (numPack - 1));
                tgtInd = srcInd + i - 1;
                self.comms(:, tgtInd) = self.comms(:, srcInd);
                self.bigMat(tgtInd, :) = self.bigMat(srcInd, :);
                self.bigMat(:, tgtInd) = self.bigMat(:, srcInd);
            end
        end
        
        function UpdateComms(self)
            numPack = length(self.packs);
            for indPack = 1:(numPack - 1)
                indComm = indPack * numPack;
                self.comms(:, indComm) = self.CommBetween(indPack, numPack);
            end
            
            for indComm = ((numPack - 1) * numPack + 1):numPack^2
                indPack = indComm - (numPack - 1) * numPack;
                self.comms(:, indComm) = self.CommBetween(numPack, indPack);
            end
        end
        
        function UpdateBigMat(self)
            numPack = length(self.packs);
            updInd = [numPack * (1:(numPack - 1)), ...
                      ((numPack - 1) * numPack + 1):numPack^2];
            allInd = 1:numPack^2;
            comm = self.comms;
            self.bigMat(updInd, allInd) = comm(:, updInd)' * comm(:, allInd);
            self.bigMat(allInd, updInd) = self.bigMat(updInd, allInd)';
        end
        
        function [ok, iniCoeffUse] = InitialCoeffUse(~, tensorUse)
            numUse = size(tensorUse, 1);
            onesVec = ones(numUse, 1);
            hess = zeros(numUse, numUse);
            for i = 1:numUse
                hess(:, i) = diag(tensorUse(:, :, i, i));
            end
            hess = (hess + hess') ./ 2;
            hessL = [hess, onesVec; onesVec', 0];
            [iniCoeffUse, ~] = linsolve(hessL, [zeros(numUse, 1); 1]);
            iniCoeffUse = iniCoeffUse(1:(end - 1));
            ok = ~isnan(sum(iniCoeffUse));
        end
        
        function [grad, hess] = GradHess(~, tensorGrad, tensorHess, coeffUse)
            numUse = length(coeffUse);
            
            grad = reshape(tensorGrad, [], numUse) * coeffUse;
            grad = reshape(grad, [], numUse) * coeffUse;
            grad = reshape(grad, [], numUse) * coeffUse;
            
            hess = reshape(tensorHess, [], numUse) * coeffUse;
            hess = reshape(hess, [], numUse) * coeffUse;
            hess = reshape(hess, [], numUse);
        end
        
        function [ok, coeffUse] = NewtonSolver(self, tGrad, tHess, coeffUse)
            lagMult = 0;
            if self.maxIterNewton == 0
                ok = true;
                return;
            end
            for iter = 1:self.maxIterNewton
                [grad, hess] = self.GradHess(tGrad, tHess, coeffUse);
                hess = (hess + hess') ./ 2;
                gradL = [grad + lagMult; 0];
                onesVec = ones(length(coeffUse), 1);
                hessL = [hess, onesVec; onesVec', 0];
                [step, ~] = linsolve(hessL, gradL);
                if isnan(sum(step))
                    disp('Inversion failed')
                    ok = false;
                    return;
                end
                coeffUse = coeffUse - step(1:end-1);
                lagMult = lagMult - step(end);
                if sqrt(mean(gradL.^2)) < self.gradNormThres
                    ok = true;
                    return;
                end
            end
            disp('Newton did not converge');
            ok = false;
        end
        
    end
    
end
