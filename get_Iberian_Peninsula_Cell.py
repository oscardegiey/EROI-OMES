#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Class_Producers import *
from Class_Storage import * 
from Class_Cell import Cell
import numpy as np

# year is a string like "2015" or "2016,...
def get_Iberian_Peninsula_Cell(re_share, import_share, nuc_lf, size, dt, year):
     
    """ Installed capacities in W"""
    
    CCGT = 37834e6
    Coal = 11291e6
    Nuclear = 7117e6
    Solar = 7187e6
    WT_Onshore = 27907e6
    WT_Offshore = 5e6
    Hydro_Dam = 18633.8e6
    Hydro_River = 5562.5e6
    PHES_Turb = 7.637e9
    PHES_Pump = 7.138e9
    
    
    """Hydro Inflow"""
    annual_prod_dam = 14794.00e9 #Wh
    annual_prod_runriver = 9107.00e9 #Wh
    
    """ Hydro Storage installed capacity in Wh"""
    PHES_Storage = 110.77e9
    Dam_Storage = annual_prod_dam/2
    
        
    """Load Factor and potential WT"""
    Load_factor_onshore_Spain =  0.2377
    Load_factor_onshore_Portugal = 0.2685
    Pot_onshore_Spain = 242965.72
    Pot_onshore_Portugal = 41071.27
    
    Pot_onshore_Cell = Pot_onshore_Spain+Pot_onshore_Portugal
    fac_WT_onshore = 4 #MW/km2
    Lim_WT_onshore = Pot_onshore_Cell*fac_WT_onshore #MW
    
    
    Load_factor_offshore_Spain =  0.38 # Mean between Atlantic and mediterranean ninja
    Load_factor_offshore_Portugal = 0.4 # TODO: Etre sur de ces valeurs
    Pot_offshore_Spain = 15564.343
    Pot_offshore_Portugal = 4361.192
    
    Pot_offshore_Cell = Pot_offshore_Spain+Pot_offshore_Portugal
    fac_WT_offshore = 4 #MW/km2
    Lim_WT_offshore = Pot_offshore_Cell*fac_WT_offshore #MW
        
    """PV potential"""
    Load_factor_pv_Spain = 0.16761826
    Load_factor_pv_Portugal = 0.16597264
    Pot_pv_Spain = 28082.67e6
    Pot_pv_Portugal = 4744.00e6
    
    Surf_disp_sol = Pot_pv_Spain + Pot_pv_Portugal 
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
    irr_file = 'MalagaGHI_Irr.txt' #TODO: To change with GHI
    wind_file = 'Faro_wind.txt'
    


    
    
    # Creer cellule
    cell = Cell("Iberian_Peninsula", re_share, import_share, nuc_lf, size,dt) #(name, re_share, imp_share, size, dt)
    
    # Set consumption
    
    cell.set_cons(year)
    
    average_cons = sum(cell.cons)/size
    
    # Set producers
    
    """ PV PARK """
    
    load_factor_pv = (Load_factor_pv_Spain*Pot_pv_Spain + \
    Load_factor_pv_Portugal*Pot_pv_Portugal)/Surf_disp_sol

    pv_park = PVPark('Iberian_Peninsula', size, dt) #(name, size, dt)
    pv_park.set_spec(Solar,Surf_lim_sol,45,Cost_PV,30,0.1,load_factor_pv)  
    
    # (Pmin,maxM2,CO2,cost,life,eff) 
    pv_park.set_irr(year,cell.name,irr_file)
    pv_park.set_prod()
 
    
    """ WT Park """
    
    Load_factor_onshore = (Load_factor_onshore_Spain*Pot_onshore_Spain+Load_factor_onshore_Portugal*Pot_onshore_Portugal)/Pot_onshore_Cell
    
    wt_park_onshore = WTPark('Iberian_Peninsula_onshore', size, dt)
    wt_park_onshore.set_spec(WT_Onshore/5e6,Lim_WT_onshore/5,12,Cost_WT_onshore,25,3.5, 25, 14, 5e6, Load_factor_onshore) 
    #  (n_min,n_max,CO2,cost,life, cut_in, cut_off, wind_maxpower, p_nom, load)
    wt_park_onshore.set_wind(year,cell.name,wind_file)
    wt_park_onshore.set_prod()
    
    Load_factor_offshore = (Load_factor_offshore_Spain*Pot_offshore_Spain+Load_factor_offshore_Portugal*Pot_offshore_Portugal)/Pot_offshore_Cell
    
    wt_park_offshore = WTPark('Iberian_Peninsula_offshore', size, dt)
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