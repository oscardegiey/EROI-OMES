#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 18:25:19 2018

@author: thibaultmartinelle
"""

import numpy as np
import scipy
from scipy.optimize import fsolve
from scipy.interpolate import CubicSpline
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

class PVPark:
    def __init__(self, region, specs, size, dt):
        self.name = region['name']
        self.size = size
        self.dt   = dt

        # Specs
        self.Pmin       = region['installed_cap']
        self.efficiency = specs['efficiency']
        self.minM2      = self.Pmin / (1000 * self.efficiency)
        self.maxM2      = region['potential'] / specs['pps']
        self.CO2perMWh  = specs['CO2perMWh']
        self.costpermsq = specs['cost']
        self.lifetime   = specs['lifetime']*8760 # hour
        self.load_fac   = region['load_factor']

        # Set irradiance
        self.irr = np.loadtxt(open(region['data_file'],'r'), skiprows=0)
        self.irr = self.irr.reshape(-1,self.dt).mean(axis=1) # Reshape according to the time step

        # Set production per m^2
        Pnom = 1000 # W/m^2, Nominal power of PV in Standard Test Conditions

        # Find the coefficient to apply to the irradiance to match the load factor
        correction_coeff = fsolve(lambda coef: self.load_fac - sum(self.irr*coef*self.efficiency)/(self.efficiency*Pnom*self.size),1)
        self.irr = self.irr * correction_coeff
        self.prod_unit = self.irr*self.efficiency

    def set_opti(self, surf_pv, p_gen):
        # Set the optimal values
        self.surf_pv = surf_pv
        self.p_gen   = p_gen

class WTPark:
    def __init__(self, region, specs, size, dt, onshore):
        if onshore:
            self.onshore = True; self.offshore = False
            self.name = region['name'] + '_onshore'
        else:
            self.onshore = False; self.offshore = True
            self.name = region['name'] + '_offshore'
        self.size = size
        self.dt   = dt

        # Specs
        self.n_max         = region['potential'] * specs['pps'] / specs['p_nom']
        self.n_min         = region['installed_cap'] / specs['p_nom']
        self.CO2perMWh     = specs['CO2perMWh']
        self.costperWTinst = specs['cost']
        self.lifetime      = specs['lifetime']*8760 # hour
        self.cut_in        = specs['cut_in']
        self.cut_off       = specs['cut_off']
        self.rated_speed   = specs['rated_speed']
        self.wind_maxpower = self.rated_speed
        self.p_nom         = specs['p_nom']
        self.load_fac      = region['load_factor']
        self.wtpc          = np.array(specs['wtpc'])

        # Set wind
        self.wind = np.loadtxt(open(region['data_file'],'r'), skiprows=0)
        self.wind = self.wind.reshape(-1,self.dt).mean(axis=1) # Reshape according to the time step

        # Set production - Interpolates wtpc with cubic splines
        splines = CubicSpline(np.linspace(self.cut_in,self.rated_speed,len(self.wtpc)),self.wtpc)
        eta = lambda wind: ((wind >= self.rated_speed) & (wind <= self.cut_off)) \
            + ((wind >= self.cut_in) & (wind <= self.rated_speed)) * splines(wind)

        # Plot WTPC
        #w = np.linspace(0,30,100)
        #plt.plot(w,eta(w))
        #plt.show()

        correction_coeff = fsolve(lambda coef: self.load_fac - sum(self.p_nom*eta(self.wind*coef))/(self.p_nom*self.size),1)
        self.prod_unit = self.p_nom * eta(self.wind * correction_coeff)

    def set_opti(self, n_wt, p_gen):
        # Set the optimal values
        self.n_wt = n_wt
        self.p_gen = p_gen

class CCGTPlant:
    def __init__(self, specs):
        self.Pmax       = specs['p_max']
        self.Pmin       = specs['p_min']
        self.CO2perMWh  = specs['CO2perMWh']
        self.efficiency = specs['efficiency']

    def set_opti(self, p_gen, gas_import):
        self.p_gen = p_gen
        self.p_NRE = gas_import * self.efficiency
        self.p_RE = self.p_gen - self.p_NRE

class CoalPlant:
    def __init__(self, specs):
        self.Pmax       = specs['p_max']
        self.Pmin       = specs['p_min']
        self.CO2perMWh  = specs['CO2perMWh']
        self.efficiency = specs['efficiency']

    def set_opti(self, p_gen):
        self.p_gen = p_gen

class NuclearPlant:
    def __init__(self, specs, size):
        self.size      = size
        self.p_inst    = specs['installed_cap']
        self.CO2perMWh = specs['CO2perMWh']
        self.load_fac  = specs['load_factor']
        self.p_ut      = self.load_fac * self.p_inst
        self.p_gen     = self.p_ut * np.ones(size)

    def set_opti(self, p_gen):
        self.p_gen = np.linspace(p_gen, p_gen, self.size)

class DamPlant:
    def __init__(self, specs, size, dt):
        self.size           = size
        self.dt             = dt
        self.p_inst         = specs['installed_cap']
        self.capacity       = specs['capacity']
        self.annual_prod    = specs['annual_prod']
        self.CO2perMWh      = specs['CO2perMWh']
        self.efficiency_out = specs['efficiency_out']

        # Local inflow m^3/s
        self.inflow = np.loadtxt(open(specs['data_file'],'r'), skiprows=1)
        if np.sum(self.inflow) != 0:
            self.inflow = self.inflow/sum(self.inflow)*self.annual_prod
        self.local_in = self.inflow.reshape(-1,self.dt).mean(axis=1)

    def set_opti(self, p_out, energy):
        self.p_out=p_out
        self.energy = energy

class RunRiver:
    def __init__(self, specs, size, dt):
        self.size        = size
        self.dt          = dt
        self.p_inst      = specs['installed_cap']
        self.annual_prod = specs['annual_prod']
        self.CO2perMWh   = specs['CO2perMWh']
        self.efficiency  = specs['efficiency']

        # Local inflow m^3/s
        self.inflow = np.loadtxt(open(specs['data_file'],'r'), skiprows=1)
        if np.sum(self.inflow) != 0:
            self.inflow = self.inflow/sum(self.inflow)*self.annual_prod
        self.local_in = self.inflow.reshape(-1,self.dt).mean(axis=1)

    def set_opti(self,p_gen):
        self.p_gen=p_gen
