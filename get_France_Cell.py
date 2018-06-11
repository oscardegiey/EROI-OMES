#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 09:48:42 2018

@author: thibaultmartinelle
"""

from Class_Producers import *
from Class_Storage import * 
from Class_Cell import Cell
import numpy as np

# year is a string like "2015" or "2016,...
def get_France_Cell(re_share, import_share, nuc_lf, size, dt, year):
     
    """ Installed capacities in W"""
    
    CCGT = 17898e6
    Coal = 2997e6
    Nuclear = 63130e6
    Solar = 7646e6
    WT_Onshore = 13539e6
    WT_Offshore = 10e6
    Hydro_Dam = 8277e6
    Hydro_River = 10327e6
    PHES_Turb = 5.512e9
    PHES_Pump = 4.317e9

    """Hydro Inflow"""
    annual_prod_dam = 21943e9 #Wh
    annual_prod_runriver = 26655e9 #Wh
    
    """ Hydro Storage installed capacity in Wh"""
    PHES_Storage =  83.37e9
    Dam_Storage = annual_prod_dam/2 
    
        
    """Renewable potentials"""
    #TODO: all randow values
    surf_lim_wt_onshore = 276096.4642 #km2 
    fac_wt_onshore = 4 #MW/km2
    lim_wt_onshore = surf_lim_wt_onshore*fac_wt_onshore #MW
    surf_lim_wt_offshore = 61446.40829 #km2
    fac_wt_offshore = 4 #MW/km2
    lim_wt_offshore = surf_lim_wt_offshore*fac_wt_offshore/2 #MW
    Surf_disp_sol = 32658.08904 * 1e6
    fac_sol = 5
    Surf_lim_sol = Surf_disp_sol/fac_sol
    
    """ Upper storage limit in Wh"""
    #TODO: random values
    Lim_PtG = 1000e8*1e8 #W
    Lim_Batt = 10e12 #Wh
    Lim_Gas = 100e12 #Wh
    
    
    """Cost"""
    #TODO: old values
    Cost_PV = 428e3
    Cost_WT_onshore = 18713340740.74074
    Cost_WT_offshore = 28668708/3600*5e6
    Cost_Batt = 132*1.01
    ESOI_PHES = 704
    ESOI_Dam = 0 # Not used
    Cost_PtG = 100
    Cost_gas = 0.001
    #Ajouter un cout au Dam pour quand on utilise le stockage comme pour le PHES??
    
    """Load Factor"""
    Load_factor_onshore = 0.2022
    Load_factor_offshore = 0.43
    Load_factor_pv = 0.1397
    
    """Files"""
    irr_file = 'Le_HavreGHI_Irr.txt'
    wind_file = 'Saint_Brieuc_wind.txt'
    

    
    
    
    
    # Creer cellule
    cell = Cell("France", re_share, import_share, nuc_lf, size,dt) #(name, re_share, imp_share, size, dt)
    
    # Set consumption
    
    cell.set_cons(year)
    
    average_cons = sum(cell.cons)/size
    
    # Set producers
    
    """ PV PARK """

    pv_park = PVPark('France', size, dt) #(name, size, dt)
    pv_park.set_spec(Solar,Surf_lim_sol,45,Cost_PV,30,0.1,Load_factor_pv)  
    
    # (Pmin,maxM2,CO2,cost,life,eff,load_fac) 
    pv_park.set_irr(year,cell.name,irr_file)
    pv_park.set_prod()
 
    
    """ WT Park """
    
    wt_park_onshore = WTPark('France_onshore', size, dt)
    wt_park_onshore.set_spec(WT_Onshore/5e6,lim_wt_onshore/5,12,Cost_WT_onshore,25,3.5, 25, 14, 5e6, Load_factor_onshore) 
    #  (n_min,n_max,CO2,cost,life, cut_in, cut_off, wind_maxpower, p_nom, load)
    #TODO: to change
    wt_park_onshore.set_wind(year,cell.name,wind_file)
    wt_park_onshore.set_prod()
    
    wt_park_offshore = WTPark('France_offshore', size, dt)
    wt_park_offshore.set_spec(WT_Offshore/5e6,lim_wt_offshore/5,12,Cost_WT_offshore,25,3.5, 25, 14, 5e6, Load_factor_offshore) 
    #  (n_min,n_max,CO2,cost,life, cut_in, cut_off, wind_maxpower, p_nom, load)
    #TODO: to change
    wt_park_offshore.set_wind(year,cell.name,wind_file)
    wt_park_offshore.set_prod()
    
    """ CCGT Plant """
    ccgt = CCGTPlant(size)
    ccgt.set_spec(0,CCGT,490,0.59)#TODO: to change
    #ccgt.set_spec(Pmin,Pmax,CO2,eff)
    
    """ Coal Plant """
    coal = CoalPlant(size)
    coal.set_spec(0, Coal, 0, 1) #Pmin,Pmax,CO2,eff)
    #TODO: CO2?
    
    """ Nuclear Plant """
    nuc = NuclearPlant(size)
    nuc.set_spec(Nuclear,0,nuc_lf) #(p_inst,CO2, load_fac)
    #TODO: CO2?
    
    """ Dam Plant """ #TODO quelle valeur pour local in?
    dam = DamPlant(size,dt)
    dam.set_spec(Hydro_Dam, Dam_Storage, 0, 1) #set_spec(p_inst, capacity, CO2,local_in, eff_out)
    dam.set_local_in(cell.name,annual_prod_dam)

    
    """ RunRiver Plant """
    runriver = RunRiver(size,dt)
    runriver.set_spec(Hydro_River, 0, 1) #(p_inst, CO2, local_in, eff)
    runriver.set_local_in(cell.name,annual_prod_runriver)
    
    cell.set_producers([pv_park, wt_park_onshore, wt_park_offshore, ccgt, nuc, dam, runriver, coal])
    
    # Set storage
    """ BATTERIES """
    batt = Battery("Lit", size, dt) #(type, size, dt)

    batt.set_spec(0, Lim_Batt, 0.9, 1,1/5, 1/5, 4.1667e-5, 0.8, 6000, 15, Cost_Batt) 

    #(Esizemin,Esizemax, efficiency_in, efficiency_out, alpha_in, alpha_out, leak, DoD, n_cycle, lifetime, eps_inst)
  
    
    """ PHES """
    
    phes = PHES(size)

    phes.set_spec(PHES_Pump, PHES_Turb, PHES_Storage, 0,np.zeros(size), 0.85, 0.9, ESOI_PHES)

    #phes.set_spec(p_inst_pump, p_inst_out, capacity, CO2,local_in, eff_out, eff_up, ESOI)
    
    
    """ PtG """
    ptg = PtG(size, dt)
    ptg.set_spec(0, Lim_PtG, 0.63, 10000, Cost_PtG,0.15, 15)

    
    """ Gas """
    gas = Gas(size)
    gas.set_spec(0, Lim_Gas, Cost_gas, 25) 
    #(Esizemin,Esizemax, eps, life)
    
    cell.set_storages([batt, phes, ptg, gas])
    
    
    # Return cell
    
    return cell