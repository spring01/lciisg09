classdef CombineDIIS < handle
    
    properties (Access = public)
        
        init;
        final;
        
    end
    
    properties (Access = private)
        
        energies = 0;
        switchThres = 1e-2;
        
    end
    
    methods
        
        function newFockList = NewFock(self, fockList, densList, energy)
            newFockListInit = self.init.NewFock(fockList, densList, energy);
            newFockListFinal = self.final.NewFock(fockList, densList, energy);
            if abs(energy - self.energies(end)) > self.switchThres
                newFockList = newFockListInit;
            else
                newFockList = newFockListFinal;
            end
            self.energies = [self.energies, energy];
        end
        
    end
    
end
