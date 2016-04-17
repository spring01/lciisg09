classdef LISTb < handle
    
    properties (Access = private)
        
        fockInVector;
        fockOutVectors;
        deltaFockVectors;
        
        densOutVectors;
        energies;
        
    end
    
    methods
        
        function obj = LISTd(initialFockVec, numVectors, type)
            if(nargin < 2)
                numVectors = 5;
            end
            if(nargin < 3)
                type = 'r';
            end
            lenVec = numel(initialFockVec);
            if(strcmpi(type, 'r'))
                obj.fockOutVectors{1} = zeros(lenVec, numVectors);
                obj.deltaFockVectors{1} = zeros(lenVec, numVectors);
                obj.densOutVectors{1} = zeros(lenVec, numVectors);
            elseif(strcmpi(type, 'u'))
                throw(MException('LISTd:LISTd', 'not implemented yet'));
            end
            obj.energies = zeros(numVectors, 1);
            obj.fockInVector = initialFockVec;
        end
        
        function Push(obj, newFockOutVector, newDensOutVector, energy)
            % push new Fock and density in
            for spin = 1:length(obj.fockOutVectors)
                obj.fockOutVectors{spin}(:, 1:end-1) = obj.fockOutVectors{spin}(:, 2:end);
                obj.fockOutVectors{spin}(:, end) = newFockOutVector(:, spin);
                
                obj.densOutVectors{spin}(:, 1:end-1) = obj.densOutVectors{spin}(:, 2:end);
                obj.densOutVectors{spin}(:, end) = 2*newDensOutVector(:, spin);
                
                obj.deltaFockVectors{spin}(:, 1:end-1) = obj.deltaFockVectors{spin}(:, 2:end);
                obj.deltaFockVectors{spin}(:, end) = (newFockOutVector - obj.fockInVector) / 2;
            end
            obj.energies(1:end-1) = obj.energies(2:end);
            obj.energies(end) = energy;
        end
        
        function PushFockIn(obj, newFockInVector)
            obj.fockInVector = newFockInVector;
        end
        
        function [optFockVector, coeffs, useFockVectors] = OptFockVector(obj)
            useFockVectors = cell(1, length(obj.fockOutVectors));
            optFockVector = zeros(size(obj.fockOutVectors{1}, 1), length(obj.fockOutVectors));
            numVectors = sum(sum(obj.fockOutVectors{1}.^2) ~= 0);
            if(numVectors == 0 || numVectors == 1)
                for spin = 1:length(obj.fockOutVectors)
                    optFockVector(:, spin) = obj.fockOutVectors{spin}(:, end);
                    useFockVectors{spin} = obj.fockOutVectors{spin}(:, end);
                end
                if(length(useFockVectors) == 1)
                    useFockVectors = useFockVectors{1};
                end
                coeffs = 1;
                return;
            end
            
            useDeltaFockVectors = cell(1, length(obj.fockOutVectors));
            useDensVectors = cell(1, length(obj.fockOutVectors));
            for spin = 1:length(obj.fockOutVectors)
                useDeltaFockVectors{spin} = obj.deltaFockVectors{spin}(:, end-numVectors+1:end);
                useFockVectors{spin} = obj.fockOutVectors{spin}(:, end-numVectors+1:end);
                useDensVectors{spin} = obj.densOutVectors{spin}(:, end-numVectors+1:end);
            end
            useEnergies = obj.energies(end-numVectors+1:end);
            
            hessian = zeros(numVectors, numVectors);
            for spin = 1:length(obj.fockOutVectors)
                for i = 1:numVectors
                    for j = 1:numVectors
                        hessian(i, j) = hessian(i, j) + useEnergies(i) ...
                            + useDeltaFockVectors{spin}(:, i)' * (useDensVectors{spin}(:, j) - useDensVectors{spin}(:, i));
                    end
                end
            end
            hessian = hessian';
            
            onesVec = -ones(numVectors, 1);
            hessian = [ ...
                hessian, onesVec; ...
                onesVec', 0];
            diisCoefficients = hessian \ [zeros(numVectors,1); -1];
            coeffs = diisCoefficients(1:end-1);
            
            for spin = 1:length(obj.fockOutVectors)
                optFockVector(:, spin) = useFockVectors{spin} * coeffs;
            end
            
            if(length(useFockVectors) == 1)
                useFockVectors = useFockVectors{1};
            end
            
%             disp(coeffs');
        end
        
    end
    
end
