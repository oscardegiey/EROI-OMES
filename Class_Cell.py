
import numpy as np
from Class_Producers import *
from Class_Storage import *
import sys


class Cell:

    def __init__(self, name, re_share, imp_share, nuc_lf, size, dt):
        """Inits Cell"""
        self.name = name
        self.re_share = re_share
        self.imp_share = imp_share
        self.nuc_lf = nuc_lf
        self.size = size
        self.dt = dt
        self.cons = np.linspace(0.,0.,self.size)
        self.pv_parks = []
        self.wt_parks = []
        self.bat_parks = []
      
    def set_cons(self, year):
        """Set consumption.
        
        Reads a file of the form 'cellname_year_cons.txt' located in the 
        subfolder 'Input_europe/'. Set the attribute cons in function of the current
        self.dt
        """
        path = "Input_europe/"+ self.name + '/' + self.name + '_Cons/' + year + '/' + self.name + "_" + year + "_cons.txt"
        try:
            file = open(path,"r")
        except FileNotFoundError:
            print("FileNotFoundError: wrong entered consumption file name," +path+ "does not exist")
            sys.exit(1)
            
        entete = file.readline()
        i=0
        tmp = np.linspace(0,0,self.size*self.dt)
        while i<self.size*self.dt:
            try:
                tmp[i]=file.readline()
            except ValueError:
                print("End of file INPUT/"+ self.name + "_cons.txt reached \n" + 
                "-> range of iteration indices:", int(self.size * self.dt), "\n-> document size:",i)
                sys.exit(1)
            i += 1
        j = 0
        while j<self.size:
            self.cons[j] = sum(tmp[self.dt*j:self.dt*j+self.dt])/self.dt
            j = j+ 1
        file.close()
        
        self.cons = self.cons*1e6 # cons is given in MW in files
     
     
    def set_producers(self, producers): 
        """Set the Cell producers.
        
        Take a list of objects producer in argument and store this object in the 
        structure of the Cell."""
        for producer in producers:
            if type(producer) is PVPark:
                self.pv_parks.append(producer)
            
            elif type(producer) is WTPark:
                self.wt_parks.append(producer)
                
            elif type(producer) is CCGTPlant:
                self.ccgt = producer
                
            elif type(producer) is NuclearPlant:
                self.nuclear = producer
                
            elif type(producer) is CoalPlant:
                self.coal = producer
            
            elif type(producer) is DamPlant:
                self.dam = producer
        
            elif type(producer) is RunRiver:
                self.runriver = producer
                
            else:
                print("Unexpected Producer type")

           
    def set_storages(self, storages):
        """Set Cell storages.
            
        Take an object storage in argument and store this object in the 
        structure of the Cell. 
        """
    
        for storage in storages:
            if type(storage) is Battery: 
                self.bat_parks.append(storage)
            
            elif type(storage) is PHES:
                self.phes = storage
            
            elif type(storage) is PtG:
                self.ptg = storage
                
            elif type(storage) is Gas:
                self.gas = storage
            
            else: 
                print("Unexpected Storage type")
            
        