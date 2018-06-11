#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 18:25:19 2018

@author: thibaultmartinelle
"""

import numpy as np
import sys


""" PV PARK """

class PVPark:
    def __init__(self, name, size, dt):
        self.name = name
        self.size = size
        self.dt = dt
        self.prod_unit = np.linspace(0.,0.,size)
        self.irr = np.linspace(0.,0.,self.size*self.dt)
    
    def set_spec(self,Pmin,maxM2,CO2,cost,life,eff,load_fac):    
        self.Pmin = Pmin #Pmin in Wp
        self.minM2 = self.Pmin/(1000*eff)
        self.maxM2 = maxM2
        self.CO2perMWh    = CO2
        self.costpermsq = cost
        self.lifetime     = life*8760
        self.efficiency   = eff
        self.load_fac = load_fac
    
    def set_irr(self,year,cellname,filename):
        """Set the PVPark irradiance.
        
        Read a file of of name filename located in subfolder
        'Input_europe/' and set the park irradiance.
        """
        path = "Input_europe/"+ cellname + '/' + cellname + '_Irr/TXT/'+ filename
        try:
            doc = open(path, "r") 
        except FileNotFoundError:
            print("FileNotFoundError: wrong entered irradiance file name, "+path+" does not exist")
            sys.exit(1)
        #entete = doc.readline()
        i=0
        while i<self.size*self.dt:
            try:
                self.irr[i]=doc.readline()
            except ValueError:
                print("End of file "+path+" reached \n" + 
                "-> range of iteration indices:", int(self.size * self.dt), "\n-> document size:",i)
                sys.exit(1)
            i = i+1
        doc.close()
    

    def set_prod(self):
        """Set the production of 1 square meter of PV.
        
        Read the irradiance and set the production according to the efficiency
        of the panel and the self.dt
        """
        
        err = 2e-3
        prod_unit_full = np.zeros(self.size*self.dt)
        while np.abs(err)>1e-3:

            self.irr=self.irr*(1+err)
            prod_unit_full = self.efficiency*self.irr
            load_tmp = sum(prod_unit_full)/(self.efficiency*1000*self.size*self.dt) # p_nom=1000W/m^2*efficiency
            err = self.load_fac - load_tmp
        
        i = 0
        while i < self.size:
            irr_av = sum(self.irr[self.dt*i:self.dt*i+self.dt])/self.dt
            self.prod_unit[i] = irr_av*self.efficiency
            i = i+1
        
    def set_opti(self, surf_pv, p_gen):
        """Set the optimal values"""
        self.surf_pv = surf_pv
        self.p_gen = p_gen
        
""" WT Park """
        
class WTPark:
    def __init__(self, name, size, dt):
        self.name = name
        self.size = size
        self.dt = dt   
        self.prod_unit = np.linspace(0.,0.,size)
        self.wind = np.linspace(0.,0.,self.size*self.dt)
        self.LF = np.linspace(0.,0.,self.size*self.dt)
    
    def set_spec(self,n_min,n_max,CO2,cost,life, cut_in, cut_off, wind_maxpower, p_nom, load_fac):    
        self.n_max = n_max
        self.n_min = n_min
        self.CO2perMWh    = CO2 #kg_CO2_eq per MWh produced
        self.costperWTinst = cost
        self.lifetime     = life*8760
        self.cut_in = cut_in
        self.cut_off = cut_off
        self.wind_maxpower = wind_maxpower
        self.p_nom = p_nom
        self.load_fac = load_fac
        
    def set_wind(self,year,cellname,filename):
        """Set the WTPark wind.
        
        Read a file of name filename located in subfolder
        'Input_europe/' and set the park wind.
        """
        
        path = "Input_europe/"+ cellname + '/' + cellname + '_Wind/TXT/'+ filename
        try:
            file = open(path,"r")
        except FileNotFoundError:
            print("FileNotFoundError: wrong entered wind file name "+path+" does not exist")
            sys.exit(1)
        #entete = file.readline()
        i = 0
        while i < self.size*self.dt:
            try:
                self.wind[i] = file.readline() 
            except ValueError:
                print("End of file "+path+" reached \n" + 
                "-> range of iteration indices:", int(self.size * self.dt), "\n-> document size:",i)
                sys.exit(1)
            i = i + 1
        file.close()
    
    def set_prod(self):
        """Set the production of the specified wind turbine.
        
        Scale the wind for a given load factor, and compute the wind production
        """
        err = 2e-3
        prod_unit_full = np.zeros(self.size*self.dt)
        while err>1e-3:
            i = 0
            while i < self.size*self.dt:
                self.wind[i]=self.wind[i]*(1+err)
                eta = self.wind[i]*0+((self.wind[i] >= self.wind_maxpower).astype(int) * (self.wind[i] <= self.cut_off).astype(int)) + ((self.wind[i] < self.wind_maxpower).astype(int) * (self.wind[i] >= self.cut_in).astype(int))*(self.wind[i]**3-self.cut_in**3)/(self.wind_maxpower**3-self.cut_in**3)
                prod_unit_full[i] = self.p_nom*eta
                i = i+1
            load_tmp = sum(prod_unit_full)/(self.p_nom*self.size*self.dt)
            err = self.load_fac - load_tmp
        
        i = 0    
        while i<self.size:
            self.prod_unit[i] = sum(prod_unit_full[self.dt*i:self.dt*i+self.dt])/self.dt
            i+=1
        
    def set_prod_LF(self, p_nom):
        """Set the production of the specified wind turbine based on the load factor file.

        """
        self.p_nom = p_nom
        path = "Input_europe/"+ self.name + '/' + self.name + '_Wind/TXT/'+ self.name + "_" + year + "_wind_LF.txt"
        try:
            file = open(path,"r")
        except FileNotFoundError:
            print("FileNotFoundError: wrong entered wind file name "+path+" does not exist")
            sys.exit(1)
        entete = file.readline()
        i = 0
        while i < self.size*self.dt:
            try:
                self.LF[i] = file.readline() 
            except ValueError:
                print("End of file "+path+" reached \n" + 
                "-> range of iteration indices:", int(self.size * self.dt), "\n-> document size:",i)
                sys.exit(1)
            i = i + 1
        file.close()
        
        while i<self.size:
            self.prod_unit[i] = sum(self.LF[self.dt*i:self.dt*i+self.dt])*self.p_nom/self.dt
    
    def set_opti(self, n_wt, p_gen):
        """Set the optimal values"""
        self.n_wt = n_wt
        self.p_gen = p_gen

        
""" CCGT Plant """

class CCGTPlant:
    def __init__(self, size):
        self.size = size
        
    def set_spec(self,Pmin,Pmax,CO2,eff):
        self.Pmax=Pmax
        self.Pmin=Pmin
        self.CO2perMWh=CO2
        self.efficiency=eff
        
    def set_opti(self, p_gen):
        self.p_gen = p_gen
        
""" Coal Plant """

class CoalPlant:
    def __init__(self, size):
        self.size = size
        
    def set_spec(self,Pmin,Pmax,CO2,eff):
        self.Pmax=Pmax
        self.Pmin=Pmin
        self.CO2perMWh=CO2
        self.efficiency=eff
        
    def set_opti(self, p_gen):
        self.p_gen = p_gen
        
""" Nuclear Plant """
        
class NuclearPlant:
    def __init__(self, size):
        self.size=size
        
    def set_spec(self,p_inst,CO2, load_fac):    
        self.p_inst = p_inst
        self.CO2perMWh=CO2
        self.load_fac = load_fac
        self.p_ut = p_inst*load_fac

    def set_prod(self):
        self.p_gen = np.linspace(self.p_inst*self.load_fac, self.p_inst*self.load_fac, self.size)
        
    def set_opti(self, p_gen):
        self.p_gen = np.linspace(p_gen, p_gen, self.size)
        

""" Hydro Plant """
        
class DamPlant:
    def __init__(self, size,dt):
        self.size = size
        self.dt = dt
        self.local_in = np.zeros(size)
        
    def set_spec(self, p_inst, capacity, CO2, eff_out):    
        self.p_inst = p_inst
        self.capacity = capacity
        self.CO2perMWh=CO2
        self.efficiency_out = eff_out
        
    def set_local_in(self,cellname,annual_prod):
        """Set the DamPlant local inflow.
        
        Read a file of name filename located in subfolder
        'Input_europe/' and set the park irradiance.
        """
        path = "Input_europe/"+ cellname + '/' + cellname + '_Local_in/'+ cellname + '_Local_in.txt'
        try:
            doc = open(path, "r") 
        except FileNotFoundError:
            print("FileNotFoundError: wrong entered local inflow file name, "+path+" does not exist")
            sys.exit(1)
        entete = doc.readline()
        i=0
        local_inflow = np.zeros(self.size*self.dt)
        while i<self.size*self.dt:
            try:
                local_inflow[i]=doc.readline()
            except ValueError:
                print("End of file "+path+" reached \n" + 
                "-> range of iteration indices:", int(self.size * self.dt), "\n-> document size:",i)
                sys.exit(1)
            i = i+1
        doc.close()
        
        integral = sum(local_inflow)
        local_inflow = local_inflow/integral*annual_prod
    
        i = 0
        while i < self.size:
            self.local_in[i] = sum(local_inflow[self.dt*i:self.dt*i+self.dt])/self.dt
            i = i+1
        
    def set_opti(self, p_out, energy): 
        self.p_out=p_out 
        #self.p_in = p_in  
        self.energy = energy 
        
        
""" Dam RunRiver Plant """
        
class RunRiver:
    def __init__(self, size, dt):
        self.size = size
        self.dt = dt
        self.local_in = np.zeros(size)
        
    def set_spec(self, p_inst, CO2, eff):    
        self.p_inst = p_inst
        self.CO2perMWh=CO2
        self.efficiency = eff
        
    def set_local_in(self,cellname,annual_prod):
        """Set the Runriver local inflow.
        
        Read a file of name filename located in subfolder
        'Input_europe/' and set the park irradiance.
        """
        path = "Input_europe/"+ cellname + '/' + cellname + '_Local_in/'+ cellname + '_Local_in.txt'
        try:
            doc = open(path, "r") 
        except FileNotFoundError:
            print("FileNotFoundError: wrong entered local inflow file name, "+path+" does not exist")
            sys.exit(1)
        entete = doc.readline()
        i=0
        local_inflow = np.zeros(self.size*self.dt)
        while i<self.size*self.dt:
            try:
                local_inflow[i]=doc.readline()
            except ValueError:
                print("End of file "+path+" reached \n" + 
                "-> range of iteration indices:", int(self.size * self.dt), "\n-> document size:",i)
                sys.exit(1)
            i = i+1
        doc.close()
        
        integral = sum(local_inflow)
        local_inflow = local_inflow/integral*annual_prod
    
        i = 0
        while i < self.size:
            self.local_in[i] = sum(local_inflow[self.dt*i:self.dt*i+self.dt])/self.dt
            i = i+1
        
    def set_opti(self,p_gen): 
        self.p_gen=p_gen
        
        
        
    
    
        
        