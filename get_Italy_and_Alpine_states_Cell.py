#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Class_Producers import *
from Class_Storage import * 
from Class_Cell import Cell
import numpy as np

# year is a string like "2015" or "2016,...
def get_Italy_and_Alpine_states_Cell(re_share, import_share, nuc_lf, size, dt, year):
     
    """ Installed capacities in W"""
    
    CCGT = 55042e6
    Coal = 8653e6
    Nuclear = 3373e6
    Solar = 21325e6
    WT_Onshore = 12178e6
    WT_Offshore = 0e6
    Hydro_Dam = 25861.15e6
    Hydro_River = 22287.07e6
    PHES_Turb = 14.175e9
    PHES_Pump = 12.398e9
    
    """Inflow"""
    annual_prod_dam = 38346.00e9 #Wh
    annual_prod_runriver = 68861.00e9 #Wh

    
    """ Hydro Storage installed capacity in Wh"""
    PHES_Storage = 512.16e9
    Dam_Storage = annual_prod_dam/2
    
        
    """Load Factor and potential WT"""
    Load_factor_onshore_Italy =  0.2042
    Load_factor_onshore_Switzerland = 0.2511
    Load_factor_onshore_Austria = 0.2728
    Pot_onshore_Italy = 134825.15
    Pot_onshore_Switzerland = 16912.14
    Pot_onshore_Austria = 29123.83
    Pot_onshore_Cell = Pot_onshore_Italy+Pot_onshore_Switzerland+Pot_onshore_Austria
    fac_WT_onshore = 4 #MW/km2
    Lim_WT_onshore = Pot_onshore_Cell*fac_WT_onshore #MW
    
    
    Load_factor_offshore_Italy =  0.331 #renewable ninja
    Load_factor_offshore_Switzerland = 0
    Load_factor_offshore_Austria = 0
    Pot_offshore_Italy = 26011.033
    Pot_offshore_Switzerland = 0
    Pot_offshore_Austria = 0
    Pot_offshore_Cell = Pot_offshore_Italy+Pot_offshore_Switzerland+Pot_offshore_Austria
    fac_WT_offshore = 4 #MW/km2
    Lim_WT_offshore = Pot_offshore_Cell*fac_WT_offshore #MW
    
        
    """PV potential"""
    Load_factor_pv_Italy = 0.15494488
    Load_factor_pv_Austria = 0.13887403
    Load_factor_pv_Switzerland = 0.15453642
    Pot_pv_Italy = 15184.28e6
    Pot_pv_Austria = 2891.07e6
    Pot_pv_Switzerland = 1706.82e6
    
    Surf_disp_sol = Pot_pv_Italy + Pot_pv_Austria + Pot_pv_Switzerland  
    fac_sol = 5
    Surf_lim_sol = Surf_disp_sol/fac_sol
    
    """ Upper storage limit in Wh"""
    #TODO: random values
    Lim_PtG = 1000e8*1e8 #W
    Lim_Batt = 10e12*1e7 #Wh
    Lim_Gas = 100e12*1e8 #Wh
    
    
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
    irr_file = 'AmendolaraGHI_Irr.txt' #TODO: To change with GHI
    wind_file = 'Vienne_wind.txt'
    

    
    
    # Creer cellule
    cell = Cell("Italy_and_Alpine_states", re_share, import_share, nuc_lf, size,dt) #(name, re_share, imp_share, size, dt)
    
    # Set consumption
    
    cell.set_cons(year)
    
    average_cons = sum(cell.cons)/size
    
    # Set producers
    
    """ PV PARK """

    load_factor_pv = (Load_factor_pv_Italy*Pot_pv_Italy + \
    Load_factor_pv_Austria*Pot_pv_Austria + Load_factor_pv_Switzerland*Pot_pv_Switzerland)/Surf_disp_sol

    pv_park = PVPark('Italy_and_Alpine_states', size, dt) #(name, size, dt)
    pv_park.set_spec(Solar,Surf_lim_sol,45,Cost_PV,30,0.1,load_factor_pv)  
    
    # (Pmin,maxM2,CO2,cost,life,eff) 
    pv_park.set_irr(year,cell.name,irr_file)
    pv_park.set_prod()
 
    
    """ WT Park """
    
    Load_factor_onshore = (Load_factor_onshore_Austria*Pot_onshore_Austria+Load_factor_onshore_Italy*Pot_onshore_Italy+Load_factor_onshore_Switzerland*Pot_onshore_Switzerland)/Pot_onshore_Cell
    
    wt_park_onshore = WTPark('Italy_and_Alpine_states_onshore', size, dt)
    wt_park_onshore.set_spec(WT_Onshore/5e6,Lim_WT_onshore/5,12,Cost_WT_onshore,25,3.5, 25, 14, 5e6, Load_factor_onshore) 
    #  (n_min,n_max,CO2,cost,life, cut_in, cut_off, wind_maxpower, p_nom, load)
    wt_park_onshore.set_wind(year,cell.name,wind_file)
    wt_park_onshore.set_prod()
    
    Load_factor_offshore = (Load_factor_offshore_Austria*Pot_offshore_Austria+Load_factor_offshore_Italy*Pot_offshore_Italy+Load_factor_offshore_Switzerland*Pot_offshore_Switzerland)/Pot_offshore_Cell
    
    wt_park_offshore = WTPark('Italy_and_Alpine_states_offshore', size, dt)
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
    #ptg.set_spec(Pmin, Pmax,  efficiency_in, n_cycle, eps_inst, alpha, life)
    
    """ Gas """
    gas = Gas(size)
    gas.set_spec(0, Lim_Gas, Cost_gas, 25) 
    #(Esizemin,Esizemax, eps, life)
    
    cell.set_storages([batt, phes, ptg, gas])
    
    
    # Return cell
    
    return cell