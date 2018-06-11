#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 22 17:46:10 2018

@author: thibaultmartinelle
"""
from pulp import *

def compute_val(cells, grid):
    dt= cells[0].dt
    size = cells[0].size
    E_prod = {}
    E_re_prod = {}
    E_imp = {}
    E_nuc = {}
    E_coal = {}
    E_ccgt = {}
    E_pv = {}
    E_wt = {}
    E_dam = {}
    E_runriver = {}
    curtailment= {}
    E_cons = {}
    
    EROI_PV = {}
    EROI_WT = {}
    Share = {}
    LF = {}
    for cell in cells:
        E_imp[cell.name] = sum(grid.e_import[cell.name]*dt)
        E_nuc[cell.name] = sum(cell.nuclear.p_gen*dt)
        E_coal[cell.name] = sum(cell.coal.p_gen*dt)
        E_ccgt[cell.name] = sum(cell.ccgt.p_gen*dt)
        E_cons[cell.name] = sum(cell.cons*dt)
        E_pv[cell.name] = 0
        E_wt[cell.name] = 0
        for park in cell.pv_parks:
            E_pv[cell.name] += sum(park.p_gen*dt)
        for park in cell.wt_parks:
            E_wt[cell.name] += sum(park.p_gen*dt)
        E_dam[cell.name] = sum(cell.dam.p_out*dt)
        E_runriver[cell.name] = sum(cell.runriver.p_gen*dt)
        E_re_prod[cell.name] =E_dam[cell.name] + E_runriver[cell.name]+E_pv[cell.name]+E_wt[cell.name]
        E_prod[cell.name] = E_re_prod[cell.name]+grid.gas_import[cell.name]*cell.ccgt.efficiency*cell.size*cell.dt +E_coal[cell.name]+E_nuc[cell.name]
        curtailment[cell.name]= sum(sum(park.p_gen for park in cell.pv_parks) + sum(park.p_gen for park in cell.wt_parks) + cell.ccgt.p_gen + cell.coal.p_gen +cell.nuclear.p_gen + \
            sum(park.p_in for park in cell.bat_parks) + sum(park.p_out for park in cell.bat_parks) + cell.phes.p_in + cell.phes.p_out+ \
            cell.dam.p_out + cell.runriver.p_gen + cell.ptg.p_in - sum(grid.e_exch[cell.name][neighbour.name] for neighbour in cells) + \
            sum(grid.e_exch[neighbour.name][cell.name] for neighbour in cells) + grid.e_import[cell.name]  - cell.cons)*dt
        
        EROI_PV[cell.name] = {}
        EROI_WT[cell.name] = {}
        Share[cell.name] = {}
        LF[cell.name] = {} 
        
        
        #Compute EROI
        for park in cell.pv_parks:
            EROI_PV[cell.name][park.name] = sum(park.prod_unit*dt)/park.costpermsq*park.lifetime/(size*dt)
            
        EROI_WT[cell.name] = {}
        for park in cell.wt_parks:
            EROI_WT[cell.name][park.name] = sum(park.prod_unit)*dt/park.costperWTinst*park.lifetime/(size*dt)
        
        #Compute Shares
        Share[cell.name]['RE'] = value(E_re_prod[cell.name])/value(E_prod[cell.name])
        Share[cell.name]['RE/cons'] = value(E_re_prod[cell.name])/value(E_cons[cell.name])
        Share[cell.name]['PV'] = value(E_pv[cell.name])/value(E_prod[cell.name])
        Share[cell.name]['WT'] = value(E_wt[cell.name])/value(E_prod[cell.name])
        Share[cell.name]['Runriver'] = value(E_runriver[cell.name])/value(E_prod[cell.name])
        Share[cell.name]['Dam'] = value(E_dam[cell.name])/value(E_prod[cell.name])
        Share[cell.name]['Coal'] = value(E_coal[cell.name])/value(E_prod[cell.name])
        Share[cell.name]['Nuc'] = value(E_nuc[cell.name])/value(E_prod[cell.name])
        Share[cell.name]['CCGT'] = value(E_ccgt[cell.name])/value(E_prod[cell.name])
        Share[cell.name]['Import'] = value(E_imp[cell.name])/value(E_prod[cell.name])
        Share[cell.name]['Curt'] = value(curtailment[cell.name])/value(E_re_prod[cell.name])
        #Compute load factor
        LF[cell.name]['Runriver'] = value(E_runriver[cell.name])/(cell.runriver.p_inst*size*dt)
        LF[cell.name]['Dam'] = value(E_dam[cell.name])/((cell.dam.p_inst*size*dt) or 1) # avoid div 0
        LF[cell.name]['PtG'] = -sum(cell.ptg.p_in*dt)/(cell.ptg.p_inst*size*dt)
        LF[cell.name]['Batt'] = -sum(cell.bat_parks[0].p_in*dt)/(cell.bat_parks[0].capacity*cell.bat_parks[0].alpha_in*size*dt)
        LF[cell.name]['CCGT'] = sum(cell.ccgt.p_gen*dt)/(cell.ccgt.Pmax*size*dt)
        
    Share['tot'] = {}
    E_prod['tot'] = 0
    E_re_prod['tot'] = 0
    E_imp['tot'] = 0
    for cell in cells:
        E_prod['tot'] += value(E_prod[cell.name])
        E_imp['tot'] += E_imp[cell.name]
        E_re_prod['tot'] += value(E_re_prod[cell.name])
    Share['tot']['RE'] = E_re_prod['tot']/(E_prod['tot']+E_imp['tot'])
    
    #------- E INVESTED ----------
    E_invested = {}
    E_invested['REassets'] = {}
    E_invested['Storage'] = {}
    E_invest_tot = {}
    E_invest_tot['REassets'] = 0
    E_invest_tot['Storage'] = 0
    storage_weight = 1/2
    E_invested_0 = {}
    for cell in cells:
        dt = cell.dt
        size = cell.size
        E_invested_0[cell.name] = 0
        E_invested['REassets'][cell.name] = 0
        E_invested['Storage'][cell.name] = 0
        E_invested['Storage'][cell.name] += -sum(cell.phes.p_in)*dt/(cell.phes.ESOI)
        E_invested['Storage'][cell.name] += sum(cell.ptg.p_in)*cell.ptg.eps_ut*(1-storage_weight)
        E_invested['Storage'][cell.name] += cell.ptg.p_inst*cell.ptg.eps_inst*(size*dt)/cell.ptg.lifetime*(storage_weight) 

        
        for park in cell.pv_parks:
            E_invested['REassets'][cell.name] += park.surf_pv*park.costpermsq*(size*dt)/park.lifetime
            E_invested_0[cell.name] += park.minM2*park.costpermsq*(size*dt)/park.lifetime
            
        for park in cell.wt_parks:
            E_invested['REassets'][cell.name] += park.n_wt*park.costperWTinst*(size*dt)/park.lifetime
            E_invested_0[cell.name] += park.n_min*park.costperWTinst*(size*dt)/park.lifetime
            
        for park in cell.bat_parks:
            E_invested['Storage'][cell.name] += park.capacity*park.eps_inst*(storage_weight)*(size*dt)/park.lifetime
            E_invested['Storage'][cell.name] += sum(park.p_in)*park.eps_ut*(1-storage_weight)
    
        E_invested['Storage'][cell.name] += cell.gas.capacity*cell.gas.eps*(size*dt)/cell.gas.lifetime 
     
        
        E_invest_tot['Storage'] += E_invested['Storage'][cell.name]
        E_invest_tot['REassets'] += E_invested['REassets'][cell.name]
    
    return [EROI_PV, EROI_WT, Share,LF, E_prod, E_invested, E_invest_tot, E_invested_0]