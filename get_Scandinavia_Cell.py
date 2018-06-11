#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Class_Producers import *
from Class_Storage import * 
from Class_Cell import Cell
import numpy as np

# year is a string like "2015" or "2016,...
def get_Scandinavia_Cell(re_share, import_share, nuc_lf, size, dt, year):
     
    """ Installed capacities in W"""
    
    CCGT = 13033e6
    Coal = 4847e6
    Nuclear = 8603e6
    Solar = 743e6
    WT_Onshore = 12064e6
    WT_Offshore = 1468e6
    Hydro_Dam = 47961e6
    Hydro_River = 6.28e6
    PHES_Turb = 1.364e9
    PHES_Pump = 0.983e9
    
    
    """Hydro Inflow"""
    annual_prod_dam = 194113.00e9 #Wh
    annual_prod_runriver = 19.00e9 #Wh
    
    """ Hydro Storage installed capacity in Wh"""
    PHES_Storage = 471.51e9
    Dam_Storage = annual_prod_dam/2 
    
        
    """Load Factor and potential WT"""
    Load_factor_onshore_Norway =  0.2873
    Load_factor_onshore_Denmark = 0.2582
    Load_factor_onshore_Sweden = 0.2361
    Pot_onshore_Norway = 171132.4
    Pot_onshore_Denmark = 19530.57
    Pot_onshore_Sweden = 102654.45
    Pot_onshore_Cell = Pot_onshore_Norway+Pot_onshore_Denmark+Pot_onshore_Sweden
    fac_WT_onshore = 4 #MW/km2
    Lim_WT_onshore = Pot_onshore_Cell*fac_WT_onshore #MW
    
    
    Load_factor_offshore_Norway =  0.4677 #Denmark LF
    Load_factor_offshore_Denmark = 0.4677 
    Load_factor_offshore_Sweden = 0.4677 #Denmark LF
    Pot_offshore_Norway = 84507.09
    Pot_offshore_Denmark = 41622.80
    Pot_offshore_Sweden = 56210.83
    Pot_offshore_Cell = Pot_offshore_Norway+Pot_offshore_Denmark+Pot_offshore_Sweden
    fac_WT_offshore = 4 #MW/km2
    Lim_WT_offshore = Pot_offshore_Cell*fac_WT_offshore #MW
    
        
    """PV potential"""
    Load_factor_pv_Sweden = 0.09758276
    Load_factor_pv_Norway = 0.09784635
    Load_factor_pv_Denmark = 0.11233852
    Pot_pv_Sweden = 8034.29e6
    Pot_pv_Norway = 17740.80e6
    Pot_pv_Denmark = 2175.82e6
    
    Surf_disp_sol = Pot_pv_Sweden + Pot_pv_Norway + Pot_pv_Denmark 
    fac_sol = 5
    Surf_lim_sol = Surf_disp_sol/fac_sol
    
    """ Upper storage limit in Wh"""
    Lim_PtG = 1000e8*1e8 #W
    Lim_Batt = 10e12*1e7 #Wh
    Lim_Gas = 100e12*1e8 #Wh
    
    
    """Cost"""
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
    irr_file = 'StockholmGHI_Irr.txt' 
    wind_file = 'Halmstad_wind.txt'
    

    
    
    # Creer cellule
    cell = Cell("Scandinavia", re_share, import_share, nuc_lf, size,dt) #(name, re_share, imp_share, size, dt)
    
    # Set consumption
    
    cell.set_cons(year)
    
    average_cons = sum(cell.cons)/size
    
    # Set producers
    
    """ PV PARK """
    
    load_factor_pv = (Load_factor_pv_Sweden*Pot_pv_Sweden + \
    Load_factor_pv_Norway*Pot_pv_Norway + \
    Load_factor_pv_Denmark*Pot_pv_Denmark)/Surf_disp_sol

    pv_park = PVPark('Scandinavia', size, dt) #(name, size, dt)
    pv_park.set_spec(Solar,Surf_lim_sol,45,Cost_PV,30,0.1,load_factor_pv)  
    
    # (Pmin,maxM2,CO2,cost,life,eff) 
    pv_park.set_irr(year,cell.name,irr_file)
    pv_park.set_prod()
 
    
    """ WT Park """
    
    Load_factor_onshore = (Load_factor_onshore_Sweden*Pot_onshore_Sweden+Load_factor_onshore_Norway*Pot_onshore_Norway+Load_factor_onshore_Denmark*Pot_onshore_Denmark)/Pot_onshore_Cell
    
    wt_park_onshore = WTPark('Scandinavia_onshore', size, dt)
    wt_park_onshore.set_spec(WT_Onshore/5e6,Lim_WT_onshore/5,12,Cost_WT_onshore,25,3.5, 25, 14, 5e6, Load_factor_onshore) 
    #  (n_min,n_max,CO2,cost,life, cut_in, cut_off, wind_maxpower, p_nom, load)
    wt_park_onshore.set_wind(year,cell.name,wind_file)
    wt_park_onshore.set_prod()
    
    Load_factor_offshore = (Load_factor_offshore_Sweden*Pot_offshore_Sweden+Load_factor_offshore_Norway*Pot_offshore_Norway+Load_factor_offshore_Denmark*Pot_offshore_Denmark)/Pot_offshore_Cell
    
    wt_park_offshore = WTPark('Scandinavia_offshore', size, dt)
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