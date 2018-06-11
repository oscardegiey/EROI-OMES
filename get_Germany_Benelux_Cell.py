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
def get_Germany_Benelux_Cell(re_share, import_share, nuc_lf, size, dt, year):
     
    """ Installed capacities in W"""
    
    CCGT = 72469e6 
    Coal = 50445e6
    Nuclear = 15914e6
    Solar = 48112e6
    WT_Onshore = 55425e6
    WT_Offshore = 7206e6
    Hydro_Dam = 698e6
    Hydro_River = 3827e6
    PHES_Turb = 9.402e9
    PHES_Pump = 8.663e9
    
    """Hydro Inflow"""
    annual_prod_dam = 633e9 #Wh
    annual_prod_runriver = 19235e9 #Wh
    
    """ Hydro Storage installed capacity in Wh"""
    PHES_Storage = 49.75e9
    Dam_Storage = annual_prod_dam/2
    
        
    """Load Factor and potential WT"""
    Load_factor_onshore_Belgium =  0.1976
    Load_factor_onshore_Lux = 0.2245
    Load_factor_onshore_Neth = 0.2278
    Load_factor_onshore_Germ = 0.19745
    Pot_onshore_Belgium = 12370.76
    Pot_onshore_Lux = 1245.60
    Pot_onshore_Neth = 15457.62
    Pot_onshore_Germ = 143315.74
    Pot_onshore_Cell = Pot_onshore_Belgium+Pot_onshore_Lux+Pot_onshore_Neth+Pot_onshore_Germ
    fac_WT_onshore = 4 #MW/km2
    Lim_WT_onshore = Pot_onshore_Cell*fac_WT_onshore #MW
    
    
    Load_factor_offshore_Belgium =  0.3698
    Load_factor_offshore_Neth = 0.4317
    Load_factor_offshore_Germ = 0.3701
    Pot_offshore_Belgium = 1393.52
    Pot_offshore_Neth = 34143.21
    Pot_offshore_Germ = 24061.69
    Pot_offshore_Cell = Pot_offshore_Belgium+Pot_offshore_Neth+Pot_offshore_Germ
    fac_WT_offshore = 4 #MW/km2
    Lim_WT_offshore = Pot_offshore_Cell*fac_WT_offshore/2 #MW
    
        
    """PV potential"""
    Load_factor_pv_Belgium = 0.12322934
    Load_factor_pv_Germ = 0.12448145
    Load_factor_pv_Neth = 0.12196376
    Load_factor_pv_Lux = 0.12969072
    Pot_pv_Belgium = 1335.24e6
    Pot_pv_Germ = 15620.93e6
    Pot_pv_Neth = 1663.84e6
    Pot_pv_Lux = 127.75e6
    Surf_disp_sol = Pot_pv_Belgium + Pot_pv_Germ + Pot_pv_Neth + Pot_pv_Lux
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
    Cost_Batt = 132
    ESOI_PHES = 704
    ESOI_Dam = 0 # Not used
    Cost_PtG = 100
    Cost_gas = 0.001
    #Ajouter un cout au Dam pour quand on utilise le stockage comme pour le PHES??
    
    
    """Files"""
    irr_file = 'Francfort1GHI_Irr.txt'
    wind_file = 'Heligoland_wind.txt'
    
    

    
    
    # Creer cellule
    cell = Cell("Germany_Benelux", re_share, import_share, nuc_lf, size,dt) #(name, re_share, imp_share, size, dt)
    
    # Set consumption
    
    cell.set_cons(year)
    
    average_cons = sum(cell.cons)/size
    
    # Set producers
    
    """ PV PARK """
    
    load_factor_pv = (Load_factor_pv_Belgium*Pot_pv_Belgium + \
    Load_factor_pv_Germ*Pot_pv_Germ + Load_factor_pv_Lux*Pot_pv_Lux + \
    Load_factor_pv_Neth*Pot_pv_Neth)/Surf_disp_sol
    
    pv_park = PVPark('Germany_Benelux', size, dt) #(name, size, dt)
    pv_park.set_spec(Solar,Surf_lim_sol,45,Cost_PV,30,0.1,load_factor_pv)  
    
    # (Pmin,maxM2,CO2,cost,life,eff) 
    pv_park.set_irr(year,cell.name,irr_file)
    pv_park.set_prod()
 
    
    """ WT Park """
    
    Load_factor_onshore = (Load_factor_onshore_Belgium*Pot_onshore_Belgium+Load_factor_onshore_Lux*Pot_onshore_Lux+Load_factor_onshore_Neth*Pot_onshore_Neth+Load_factor_onshore_Germ*Pot_onshore_Germ)/Pot_onshore_Cell
    
    wt_park_onshore = WTPark('Germay_Benelux_onshore', size, dt)
    wt_park_onshore.set_spec(WT_Onshore/5e6,Lim_WT_onshore/5,12,Cost_WT_onshore,25,3.5, 25, 14, 5e6, Load_factor_onshore) 
    #  (n_min,n_max,CO2,cost,life, cut_in, cut_off, wind_maxpower, p_nom, load)
    wt_park_onshore.set_wind(year,cell.name,wind_file)
    wt_park_onshore.set_prod()
    
    Load_factor_offshore = (Load_factor_offshore_Belgium*Pot_offshore_Belgium+Load_factor_offshore_Neth*Pot_offshore_Neth+Load_factor_offshore_Germ*Pot_offshore_Germ)/Pot_offshore_Cell
    
    wt_park_offshore = WTPark('Germay_Benelux_offshore', size, dt)
    wt_park_offshore.set_spec(WT_Offshore/5e6,Lim_WT_offshore/5,12,Cost_WT_offshore,25,3.5, 25, 14, 5e6, Load_factor_offshore) 
    #  (n_min,n_max,CO2,cost,life, cut_in, cut_off, wind_maxpower, p_nom, load)
    wt_park_offshore.set_wind(year,cell.name,wind_file)
    wt_park_offshore.set_prod()
    
    """ CCGT Plant """
    ccgt = CCGTPlant(size)
    ccgt.set_spec(0,CCGT,490,0.59)
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