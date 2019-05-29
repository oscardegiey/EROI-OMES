from pulp import *
import os, sys
#import cplex
import numpy as np
import time
import statistics as stat
import json
from colorama import Style
np.set_printoptions(threshold=np.inf)

## This file builds and solves the linear problem

def size_dict(d):
    s = 0
    for k,v in d.items():
        if isinstance(v, dict)              : s += size_dict(v)
        elif isinstance(v, pulp.LpVariable) : s += 1
        elif isinstance(v, list)            : s += len(v)
    return s


def EROI_optimisation(case):
    if case.log_level >= 1:
        print("\n" + Style.BRIGHT + "=== EROI optimisation ===" + Style.NORMAL)

    start          = time.time()
    dtperday       = int(24/case.dt)
    cells          = case.cells
    re_share       = case.policies[0]
    imp_share      = case.policies[1]
    storage_weight = 0.5
    
    curtailment = {}
    if case.optiTD:
        simpleTD = [None] * case.nrTD
        TDarray = [None] * case.size
        TDtuple = np.zeros(case.nrTD, dtype=int)
        TDays = np.zeros(case.size, dtype=int)
        hTD = np.zeros(case.size, dtype=int)
        runner = 0
        TDays = case.TDays
        
        for p in range(case.size):
            if TDays[p] not in simpleTD:
                simpleTD[runner] = TDays[p]
                runner += 1                
        simpleTD.sort()
        
        for p in range(case.size):
            TDarray[p] = simpleTD.index(TDays[p])

        for t in range (0,case.nrTD):
            for p in range (case.size):
                if int(TDarray[p])==t:
                    TDtuple[t] += 1
        TDtuple = TDtuple/dtperday
        
        for t in range(365): 
            hTD[t*dtperday:t*dtperday+dtperday] = range(dtperday)
    

    
    # Optimisation Problem : Problem
    problem = LpProblem(case.name, LpMinimize) 

    # Optimisation Problem : Variables
    ## -- 1 -- Electricity and Gas import of each cell
    e_import   = {}
    gas_import = {} # expressed in W
    for cell in cells:
        e_import[cell.name]   = [None] * case.size
        gas_import[cell.name] = [None] * case.size
        for t in range(case.size):
            e_import[cell.name][t]   = LpVariable("e_import_{0}_{1}".format(cell.name, t), 0, case.e_imp[cell.name])
            gas_import[cell.name][t] = LpVariable("gas_import_{0}_{1}".format(cell.name, t), 0, case.gas_imp[cell.name])
    
        
    ## -- 2 -- Electricity and Gas exchange between cells at every time step
    e_exch   = {}
    gas_exch = {}
    for cell in cells:
        e_exch[cell.name]   = {}
        gas_exch[cell.name] = {}
        for neighbour in cells:
            if case.optiTD:
                find  = np.zeros(case.nrTD*dtperday, dtype=int)
                e_exch[cell.name][neighbour.name] = [None] * case.nrTD * dtperday
                gas_exch[cell.name][neighbour.name] = [None] * case.nrTD * dtperday
                
                for t in range(case.size):
                    h = int(hTD[t])
                    TD = int(TDarray[t])
                
                    if find[h+dtperday*TD]==0: # Enters only 24*nrTD times instead of 8760
                        find[h+dtperday*TD] = 1
                        e_exch[cell.name][neighbour.name][h+dtperday*TD] = LpVariable("e_exch{0}_{1}_{2}".format(cell.name,neighbour.name, h+dtperday*TD), -case.e_exch[cell.name][neighbour.name], case.e_exch[cell.name][neighbour.name])
                        gas_exch[cell.name][neighbour.name][h+dtperday*TD] = LpVariable("gas_exch{0}_{1}_{2}".format(cell.name,neighbour.name, h+dtperday*TD), -case.gas_exch[cell.name][neighbour.name], case.gas_exch[cell.name][neighbour.name])
            else:
                e_exch[cell.name][neighbour.name]   = [None] * case.size
                gas_exch[cell.name][neighbour.name] = [None] * case.size
                for t in range(case.size):
                    e_exch[cell.name][neighbour.name][t]   = LpVariable("e_exch_{0}_{1}_{2}".format(cell.name,neighbour.name,t), -case.e_exch[cell.name][neighbour.name], case.e_exch[cell.name][neighbour.name])
                    gas_exch[cell.name][neighbour.name][t] = LpVariable("gas_exch_{0}_{1}_{2}".format(cell.name,neighbour.name,t), -case.gas_exch[cell.name][neighbour.name], case.gas_exch[cell.name][neighbour.name])
                               
    ## -- 3 -- Cell Variables
    cell_variables = {}
    for cell in cells:
        if case.optiTD:
            find      = np.zeros(case.nrTD*dtperday, dtype = int)
        variables = {}
        
        ### -- 3.1 -- Cell Variables: PV
        for park in cell.pv_parks:
            variables["{0}_surf_pv".format(park.name)] = LpVariable("{0}_{1}_surf_pv".format(cell.name, park.name), park.minM2, park.maxM2)
        
        ### -- 3.2 -- Cell Variables: WT
        for park in cell.wt_parks: 
            variables["{0}_n_wt".format(park.name)] = LpVariable("{0}_{1}_n_wt".format(cell.name, park.name), park.n_min, park.n_max)
            
        ### -- 3.3 -- Cell Variables: Batteries
        for park in cell.bat_parks:
            variables["{0}_capacity_bat".format(park.name)] = LpVariable("{0}_{1}_capacity_bat".format(cell.name, park.name), park.Esizemin, park.Esizemax)
        
        ### -- 3.4 -- Cell Variables: Power-to-Gas 
        variables['p_inst_ptg'] = LpVariable("{0}_p_inst_ptg".format(cell.name), cell.ptg.Pmin, cell.ptg.Pmax) 
        
        ### -- 3.5 -- Cell Variables: Gas capacity
        variables['capacity_gas'] = LpVariable("{0}_capacity_gas".format(cell.name), cell.gas.Esizemin, cell.gas.Esizemax)
        
        ### -- 3.6 -- Cell Variables: Optimised at every time step
        variables['p_nuc'] = LpVariable("{0}_p_nuc".format(cell.name), 0, cell.nuclear.p_ut)


        if case.optiTD==1:
            variables['p_coal']         = [None] * case.nrTD*dtperday # Coal - Initialisation
            variables['p_ccgt']         = [None] * case.nrTD*dtperday # CCGT - Initialisation
            variables['p_in_phes']      = [None] * case.nrTD*dtperday # PHES power in - Initialisation
            variables['p_out_dam']      = [None] * case.nrTD*dtperday # Dam power out - Initialisation
            variables['p_in_dam']       = [None] * case.nrTD*dtperday # Dam power in - Initialisation
            variables['p_gen_runriver'] = [None] * case.nrTD*dtperday # Run-of-the-river - Initialisation
            variables['p_out_phes']     = [None] * case.nrTD*dtperday # PHES power out - Initialisation
            variables['p_in_ptg']       = [None] * case.nrTD*dtperday # PtG power in - Initialisation
            for park in cell.bat_parks:
                variables["{0}_p_out_bat".format(park.name)] = [None] * case.nrTD*dtperday # Battery power out - Initialisation
                variables["{0}_p_in_bat".format(park.name)]  = [None] * case.nrTD*dtperday # Battery power in - Initialisation
                
        else: 
            variables['p_coal']         = [None] * case.size # Coal - Initialisation
            variables['p_ccgt']         = [None] * case.size # CCGT - Initialisation
            variables['p_in_phes']      = [None] * case.size # PHES power in - Initialisation
            variables['p_out_dam']      = [None] * case.size # Dam power out - Initialisation
            variables['p_in_dam']       = [None] * case.size # Dam power in - Initialisation
            variables['p_gen_runriver'] = [None] * case.size # Run-of-the-river - Initialisation
            variables['p_out_phes']     = [None] * case.size # PHES power out - Initialisation
            variables['p_in_ptg']       = [None] * case.size # PtG power in - Initialisation
            for park in cell.bat_parks:
                variables["{0}_p_out_bat".format(park.name)] = [None] * case.size # Battery power out - Initialisation
                variables["{0}_p_in_bat".format(park.name)]  = [None] * case.size # Battery power in - Initialisation
                
        for park in cell.bat_parks:        
            variables["{0}_e_bat".format(park.name)] = [None] * case.size # Battery energy stored - Initialisation  
        variables['e_dam']  = [None] * case.size # Dam energy stored - Initialisation
        variables['e_phes'] = [None] * case.size # PHES energy stored - Initialisation
        variables['e_gas']  = [None] * case.size # Gas energy stored - Initialisation

        for t in range(case.size):
            
            if case.optiTD==1:
                
                h = int(hTD[t])
                TD = int(TDarray[t])
                
                if find[h+dtperday*TD]==0: # Enters only 24*nrTD times instead of 8760
                    find[h+dtperday*TD] = 1
                
                    variables['p_coal'][h+dtperday*TD] = LpVariable("{0}_p_coal_{1}".format(cell.name, h+dtperday*TD), cell.coal.Pmin, cell.coal.Pmax) # Coal
                    variables['p_in_phes'][h+dtperday*TD] = LpVariable("{0}_p_in_phes_{1}".format(cell.name, h+dtperday*TD), - cell.phes.p_inst_in, 0) # PHES power in
                    variables['p_ccgt'][h+dtperday*TD] = LpVariable("{0}_p_ccgt_{1}".format(cell.name, h+dtperday*TD), cell.ccgt.Pmin, cell.ccgt.Pmax) # CCGT
                    variables['p_out_dam'][h+dtperday*TD] = LpVariable("{0}_p_out_dam_{1}".format(cell.name, h+dtperday*TD), 0, cell.dam.p_inst) # Dam power out
                    variables['p_in_dam'][h+dtperday*TD] = LpVariable("{0}_p_in_dam_{1}".format(cell.name, h+dtperday*TD), - cell.dam.local_in[t], 0) # Dam power in
                    variables['p_in_ptg'][h+dtperday*TD] = LpVariable("{0}_p_in_ptg_{1}".format(cell.name, h+dtperday*TD), None, 0) # PtG power in
                    variables['p_gen_runriver'][h+dtperday*TD] = LpVariable("{0}_p_gen_runriver_{1}".format(cell.name, h+dtperday*TD), 0, cell.runriver.local_in[t]*cell.runriver.efficiency) # Run-of-the-river
                    variables['p_out_phes'][h+dtperday*TD] = LpVariable("{0}_p_out_phes_{1}".format(cell.name, h+dtperday*TD), 0, cell.phes.p_inst_out) # PHES power out
                    
                    for park in cell.bat_parks:
                        variables["{0}_p_in_bat".format(park.name)][h+dtperday*TD] = LpVariable("{0}_{1}_p_in_bat_{2}".format(cell.name, park.name, h+dtperday*TD), None, 0) # Battery power in
                        variables["{0}_p_out_bat".format(park.name)][h+dtperday*TD] = LpVariable("{0}_{1}_p_out_bat_{2}".format(cell.name, park.name, h+dtperday*TD), 0, None) # Battery power out
                
            else:
                
                variables['p_coal'][t] = LpVariable("{0}_p_coal_{1}".format(cell.name, t), cell.coal.Pmin, cell.coal.Pmax) # Coal
                variables['p_ccgt'][t] = LpVariable("{0}_p_ccgt_{1}".format(cell.name, t), cell.ccgt.Pmin, cell.ccgt.Pmax) # CCGT
                variables['p_in_phes'][t] = LpVariable("{0}_p_in_phes_{1}".format(cell.name, t), - cell.phes.p_inst_in, 0) # PHES power in
                variables['p_out_dam'][t] = LpVariable("{0}_p_out_dam_{1}".format(cell.name, t), 0, cell.dam.p_inst) # Dam power out
                variables['p_in_dam'][t] = LpVariable("{0}_p_in_dam_{1}".format(cell.name, t), - cell.dam.local_in[t], 0) # Dam power in
                variables['p_in_ptg'][t] = LpVariable("{0}_p_in_ptg_{1}".format(cell.name, t), None, 0) # PtG power in
                variables['p_gen_runriver'][t] = LpVariable("{0}_p_gen_runriver_{1}".format(cell.name, t), 0, cell.runriver.local_in[t]*cell.runriver.efficiency) # Run-of-the-river # Run-of-the-river
                variables['p_out_phes'][t] = LpVariable("{0}_p_out_phes_{1}".format(cell.name, t), 0, cell.phes.p_inst_out) # PHES power out
                
                for park in cell.bat_parks:
                    variables["{0}_p_in_bat".format(park.name)][t] = LpVariable("{0}_{1}_p_in_bat_{2}".format(cell.name, park.name, t), None, 0) # Battery power in
                    variables["{0}_p_out_bat".format(park.name)][t] = LpVariable("{0}_{1}_p_out_bat_{2}".format(cell.name, park.name, t), 0, None) # Battery power out
                    
 
            variables['e_dam'][t] = LpVariable("{0}_e_dam_{1}".format(cell.name, t), 0, cell.dam.capacity) # Dam energy stored
            variables['e_phes'][t] = LpVariable("{0}_e_phes_{1}".format(cell.name, t), 0, cell.phes.capacity) # PHES energy stored 
            variables['e_gas'][t] = LpVariable("{0}_e_gas_{1}".format(cell.name, t), 0, None) # Gas energy stored
            for park in cell.bat_parks:
                variables["{0}_e_bat".format(park.name)][t] = LpVariable("{0}_{1}_e_bat_{2}".format(cell.name, park.name, t), 0, None) # Battery energy stored
       
        cell_variables[cell.name] = variables 
    if case.log_level >= 1:
        nvars = size_dict(cell_variables) + size_dict(e_exch) + size_dict(gas_exch) + size_dict(e_import) + size_dict(gas_import)
        print("Number of LpVariables : \t{0}".format(nvars))
        
    # Optimisation Problem: Objective function
    E_invested = 0
    for cell in cells:
        variables = cell_variables[cell.name]
        for t in range(case.size):  
            if case.optiTD==1:
                h = int(hTD[t])
                TD = int(TDarray[t])
                E_invested += -variables['p_in_phes'][h+dtperday*TD] * case.dt / (cell.phes.ESOI)
                E_invested += variables['p_in_ptg'][h+dtperday*TD] * cell.ptg.eps_ut * (1 - storage_weight)    
                E_invested += sum(variables["{0}_p_in_bat".format(park.name)][h+dtperday*TD] * park.eps_ut * (1 - storage_weight) for park in cell.bat_parks)
            else:   
                E_invested += -variables['p_in_phes'][t] * case.dt / (cell.phes.ESOI)
                E_invested += variables['p_in_ptg'][t] * cell.ptg.eps_ut * (1 - storage_weight)    
                E_invested += sum(variables["{0}_p_in_bat".format(park.name)][t] * park.eps_ut * (1 - storage_weight) for park in cell.bat_parks)
            				
        E_invested += sum(variables["{0}_surf_pv".format(park.name)] * park.costpermsq * (case.size * case.dt) / park.lifetime for park in cell.pv_parks)
        E_invested += sum(variables["{0}_n_wt".format(park.name)] * park.costperWTinst * (case.size * case.dt) / park.lifetime for park in cell.wt_parks)
        E_invested += sum(variables["{0}_capacity_bat".format(park.name)] * park.eps_inst * storage_weight * (case.size * case.dt) / park.lifetime for park in cell.bat_parks)
        E_invested += variables['capacity_gas']*cell.gas.eps*(case.size*case.dt)/cell.gas.lifetime 
        E_invested += variables['p_inst_ptg']*cell.ptg.eps_inst*(case.size*case.dt)/cell.ptg.lifetime*(storage_weight)
            
    problem += E_invested
    
    # Optimisation problem: Constraints    
    ## -- 1 -- Minimum RE share and maximum nuclear share in every cell production
    E_PV = {}; E_WT = {}; E_runriver = {}; E_dam = {}
    E_nuc = {}; E_coal = {}; E_ccgt_imp = {}; E_imp = {}
    E_RE = {}; E_prod = {}
    for cell in cells:
        variables = cell_variables[cell.name]
        if case.optiTD:
            E_runriver[cell.name] = 0; E_coal[cell.name] = 0; E_dam[cell.name] = 0
            E_PV[cell.name] = 0; E_WT[cell.name] = 0
            for TD in range(case.nrTD):
                for h in range(0,dtperday):
                    E_runriver[cell.name] += variables["p_gen_runriver"][h+dtperday*TD] * case.dt * TDtuple[TD]
                    E_coal[cell.name] += variables["p_coal"][h+dtperday*TD] * case.dt * TDtuple[TD]
                    E_dam[cell.name] += variables["p_out_dam"][h+dtperday*TD] * case.dt * TDtuple[TD]
                    E_PV[cell.name] += sum(variables["{0}_surf_pv".format(park.name)] * park.prod_unit[h+dtperday*TD] * case.dt for park in cell.pv_parks) * TDtuple[TD]
                    E_WT[cell.name] += sum(variables["{0}_n_wt".format(park.name)] * park.prod_unit[h+dtperday*TD] * case.dt for park in cell.wt_parks) * TDtuple[TD]
        else: 
            E_runriver[cell.name] = lpSum(variables["p_gen_runriver"]) * case.dt
            E_coal[cell.name]     = lpSum(variables["p_coal"]) * case.dt
            E_dam[cell.name]      = lpSum(variables["p_out_dam"]) * case.dt
            E_PV[cell.name]       = sum(variables["{0}_surf_pv".format(park.name)] * sum(park.prod_unit) * case.dt for park in cell.pv_parks)
            E_WT[cell.name]       = sum(variables["{0}_n_wt".format(park.name)] * sum(park.prod_unit) * case.dt for park in cell.wt_parks)
            
        
        E_nuc[cell.name]      = variables["p_nuc"] * case.dt * case.size
        E_ccgt_imp[cell.name] = lpSum(gas_import[cell.name]) * cell.ccgt.efficiency * case.dt
        E_imp[cell.name]      = lpSum(e_import[cell.name]) * case.dt
        
        E_RE[cell.name]   = E_PV[cell.name] + E_WT[cell.name] + E_runriver[cell.name] + E_dam[cell.name]
        E_prod[cell.name] = E_RE[cell.name] + E_nuc[cell.name] + E_coal[cell.name] + E_ccgt_imp[cell.name] + E_imp[cell.name]
        
        problem += (E_RE[cell.name] >= E_prod[cell.name] * cell.re_share)
        problem += (E_nuc[cell.name] <= E_prod[cell.name] * cell.nuc_share)
        
    ## -- 2 -- Minimum RE share in system production
    problem += (sum(E_RE.values()) >= sum(E_prod.values()) * re_share)
    
    ## -- 3 -- Maximum electricity import share in every cell
    E_cons = {}
    for cell in cells:
        E_cons[cell.name] = sum(cell.cons) * case.dt
        problem += (E_imp[cell.name] <= cell.imp_share * E_cons[cell.name])
    
    ## -- 4 -- Maximum electricity import share in system
    problem += (sum(E_imp.values()) <= cell.imp_share * sum(E_cons.values()))
    
    ## -- 5 -- Power balance at every time step and in every cell   
        
    for cell in cells:
        variables = cell_variables[cell.name]             
        if case.optiTD:
            find  = np.zeros(case.nrTD*dtperday, dtype=int)
        
        for t in range(case.size):
            power_balance = 0
            
            if case.optiTD==1:
                h = int(hTD[t])
                TD = int(TDarray[t])
                
                if find[h+dtperday*TD]==0:
                    find[h+dtperday*TD] = 1
                    power_balance += variables['p_ccgt'][h+dtperday*TD]
                    power_balance += variables['p_coal'][h+dtperday*TD]
                    power_balance += variables['p_in_phes'][h+dtperday*TD]
                    power_balance += variables['p_out_phes'][h+dtperday*TD]
                    power_balance += variables['p_out_dam'][h+dtperday*TD]
                    power_balance += variables['p_gen_runriver'][h+dtperday*TD]
                    power_balance += variables['p_in_ptg'][h+dtperday*TD]
                    power_balance += sum(variables["{0}_surf_pv".format(park.name)] * park.prod_unit[dtperday*(int(TDays[t])-1)+h] for park in cell.pv_parks)
                    power_balance += sum(variables["{0}_n_wt".format(park.name)] * park.prod_unit[dtperday*(int(TDays[t])-1)+h] for park in cell.wt_parks)
                    power_balance += sum(variables["{0}_p_in_bat".format(park.name)][h+dtperday*TD] for park in cell.bat_parks)
                    power_balance += sum(variables["{0}_p_out_bat".format(park.name)][h+dtperday*TD] for park in cell.bat_parks)
                    power_balance += variables['p_nuc']
                    power_balance += sum(e_exch[neighbour.name][cell.name][h+dtperday*TD] for neighbour in cells)
                    power_balance += e_import[cell.name][t]
                    power_balance -= cell.cons[dtperday*(int(TDays[t])-1)+h]
                    power_balance = power_balance*TDtuple[TD]
                    problem += (power_balance >= 0)
            else:
                power_balance += variables['p_ccgt'][t]
                power_balance += variables['p_coal'][t]
                power_balance += variables['p_in_phes'][t]
                power_balance += variables['p_out_phes'][t]
                power_balance += variables['p_out_dam'][t]
                power_balance += variables['p_gen_runriver'][t]
                power_balance += variables['p_in_ptg'][t]
                power_balance += sum(variables["{0}_surf_pv".format(park.name)] * park.prod_unit[t] for park in cell.pv_parks)
                power_balance += sum(variables["{0}_n_wt".format(park.name)] * park.prod_unit[t] for park in cell.wt_parks)
                power_balance += sum(variables["{0}_p_in_bat".format(park.name)][t] for park in cell.bat_parks)
                power_balance += sum(variables["{0}_p_out_bat".format(park.name)][t] for park in cell.bat_parks)
                
                power_balance += variables['p_nuc']
                power_balance += sum(e_exch[neighbour.name][cell.name][t] for neighbour in cells)
                power_balance += e_import[cell.name][t]
                power_balance -= cell.cons[t]
                problem += (power_balance >= 0)    
    ## -- 6 -- Storage balance
    ### -- 6.1 -- Storage balance: Batteries
    for cell in cells:
        variables = cell_variables[cell.name]
        
        for t in range(case.size):
            for park in cell.bat_parks:
                if case.optiTD==1:
                    h = int(hTD[t])
                    
                    if h==0:
                        h = dtperday-1
                        TD = int(TDarray[t-1])
                        
                    else:
                        h = h - 1
                        TD = int(TDarray[t-1])
                    e_in_bat = variables["{0}_p_in_bat".format(park.name)][h+dtperday*TD] * park.efficiency_in * case.dt
                    e_out_bat = variables["{0}_p_out_bat".format(park.name)][h+dtperday*TD] * (1 / park.efficiency_out) * case.dt
                    
                else:
                    e_in_bat = variables["{0}_p_in_bat".format(park.name)][t-1] * park.efficiency_in * case.dt
                    e_out_bat = variables["{0}_p_out_bat".format(park.name)][t-1] * (1 / park.efficiency_out) * case.dt
                    
                e_bat_now = variables["{0}_e_bat".format(park.name)][t]
                e_bat_prev = variables["{0}_e_bat".format(park.name)][t-1]
                leakage_coef = pow((1 - park.leak), case.dt)
                problem += (e_bat_now == e_bat_prev * leakage_coef - e_in_bat - e_out_bat)
                
    ### -- 6.2 -- Storage balance: PHES
    for cell in cells:
        variables = cell_variables[cell.name]
        for t in range(case.size):
            if case.optiTD==1:
                h = int(hTD[t])
                
                if h==0:
                    h = dtperday-1
                    TD = int(TDarray[t-1])
                       
                else:
                    h = h - 1
                    TD = int(TDarray[t-1])
                    
                e_in_phes = variables['p_in_phes'][h+dtperday*TD] * cell.phes.efficiency_in * case.dt
                e_out_phes = variables['p_out_phes'][h+dtperday*TD] * (1 / cell.phes.efficiency_out) * case.dt
                
            else:
                e_in_phes = variables['p_in_phes'][t-1] * cell.phes.efficiency_in * case.dt
                e_out_phes = variables['p_out_phes'][t-1] * (1 / cell.phes.efficiency_out) * case.dt
                
            e_phes_now = variables['e_phes'][t]
            e_phes_prev = variables['e_phes'][t-1]
            
            
            
            problem += (e_phes_now == e_phes_prev - e_in_phes - e_out_phes)
    
    ### -- 6.3 -- Storage balance: Dam
    for cell in cells:
        variables = cell_variables[cell.name]
        for t in range(case.size):
            if case.optiTD==1:
                h = int(hTD[t])
                
                if h==0:
                    h = dtperday-1
                    TD = int(TDarray[t-1])
                       
                else:
                    h = h - 1
                    TD = int(TDarray[t-1])
                    
                e_in_dam = variables['p_in_dam'][h+dtperday*TD] * case.dt
                e_out_dam = variables['p_out_dam'][h+dtperday*TD] * case.dt
                
            else: 
                e_in_dam = variables['p_in_dam'][t-1] * case.dt
                e_out_dam = variables['p_out_dam'][t-1] * case.dt
            e_dam_now = variables['e_dam'][t]
            e_dam_prev = variables['e_dam'][t-1]

            
            problem += (e_dam_now == e_dam_prev - e_in_dam - e_out_dam)
            
    ### -- 6.4 -- Storage balance: Gas
    for cell in cells: 
        variables = cell_variables[cell.name]
        for t in range(case.size):
            if case.optiTD==1:
                h = int(hTD[t])
                
                if h==0:
                    h = dtperday-1
                    TD = int(TDarray[t-1])
                       
                else:
                    h = h - 1
                    TD = int(TDarray[t-1])
                    
                e_out_gas = variables['p_ccgt'][h+dtperday*TD] * (1 / cell.ccgt.efficiency) * case.dt
                e_in_gas = variables['p_in_ptg'][h+dtperday*TD] * cell.ptg.efficiency_in * case.dt
                e_exch_gas = sum(gas_exch[cell.name][neighbour.name][h+dtperday*TD] * case.dt for neighbour in cells)
                
            else:
                e_out_gas = variables['p_ccgt'][t-1] * (1 / cell.ccgt.efficiency) * case.dt
                e_in_gas = variables['p_in_ptg'][t-1] * cell.ptg.efficiency_in * case.dt 
                e_exch_gas = sum(gas_exch[cell.name][neighbour.name][t-1] * case.dt for neighbour in cells)
            
            e_gas_now = variables['e_gas'][t]
            e_gas_prev = variables['e_gas'][t-1]
            e_imp_gas = gas_import[cell.name][t-1] * case.dt
            
            problem += (e_gas_now == e_gas_prev - e_in_gas - e_out_gas - e_exch_gas + e_imp_gas)
            
    ## -- 7 -- Power-to-Gas minimal use in every cell
    ptg_use = {}
    for cell in cells:
        variables = cell_variables[cell.name]
        ptg_use[cell.name] = 0
        for t in range(case.size):
            if case.optiTD==1:
                h = int(hTD[t])
                TD = int(TDarray[t])
                    
                ptg_use[cell.name] += variables['p_in_ptg'][h+dtperday*TD] * case.dt
                
            else:
                ptg_use[cell.name] += variables['p_in_ptg'][t] * case.dt
                
        problem += (ptg_use[cell.name] <= - cell.ptg.alpha * variables['p_inst_ptg'] * case.size * case.dt) 
        
    ## -- 8 -- Power exchange in both directions (electricity and gas)
    for cell1 in cells:
        for cell2 in cells:
            for t in range(case.size):
                if case.optiTD:
                    h = int(hTD[t])
                    TD = int(TDarray[t])
                    problem += (gas_exch[cell1.name][cell2.name][h+dtperday*TD] == -gas_exch[cell2.name][cell1.name][h+dtperday*TD])
                    problem += (e_exch[cell1.name][cell2.name][h+dtperday*TD] == -e_exch[cell2.name][cell1.name][h+dtperday*TD])    
                else:
                    problem += (gas_exch[cell1.name][cell2.name][t] == -gas_exch[cell2.name][cell1.name][t])
                    problem += (e_exch[cell1.name][cell2.name][t] == -e_exch[cell2.name][cell1.name][t])    
                
    ## -- 9 -- Limits
    ### -- 9.1 -- Limits: Battery power in/out and energy stored
    for cell in cells:
        variables = cell_variables[cell.name]
        for park in cell.bat_parks:
            max_power_in = variables["{0}_capacity_bat".format(park.name)] * park.alpha_in
            max_power_out = variables["{0}_capacity_bat".format(park.name)] * park.alpha_out
            min_energy = variables["{0}_capacity_bat".format(park.name)] * (1 - park.DoD)
            max_energy = variables["{0}_capacity_bat".format(park.name)]
            
            for t in range(case.size):
                if case.optiTD==1:
                    h = int(hTD[t])
                    TD = int(TDarray[t])
                    power_in = variables["{0}_p_in_bat".format(park.name)][h+dtperday*TD]
                    power_out = variables["{0}_p_out_bat".format(park.name)][h+dtperday*TD]
                    
                else:
                    power_in = variables["{0}_p_in_bat".format(park.name)][t]
                    power_out = variables["{0}_p_out_bat".format(park.name)][t]
                energy = variables["{0}_e_bat".format(park.name)][t]
                
                problem += (power_in >= -max_power_in)
                problem += (power_out <= max_power_out)
                problem += (energy >= min_energy)
                problem += (energy <= max_energy)
                
    ### -- 9.2 -- Limits: Power-to-Gas power in
    for cell in cells:
        variables = cell_variables[cell.name]
        max_power_in = variables['p_inst_ptg']
        
        for t in range(case.size):
            if case.optiTD==1:
                h = int(hTD[t])
                TD = int(TDarray[t])
                power_in = variables['p_in_ptg'][h+dtperday*TD]
                
            else:
                power_in = variables['p_in_ptg'][t]  
                
            problem += (power_in >= -max_power_in)
 
    ### -- 9.3 -- Limits: Gas energy stored
    for cell in cells:
        variables = cell_variables[cell.name]
        max_energy = variables['capacity_gas']
        
        for t in range(case.size):
            energy = variables['e_gas'][t]
            problem += (energy <= max_energy)
        
    building = time.time() - start    
    if case.log_level >= 1:
        print("Building problem : \t\t{0:.3f} s".format(time.time() - start))
    if case.log_level == 2:
        print("Solver output" + Style.DIM)
    else:
        old_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
                
        
    # Optimisation problem: Solve    
    start = time.time()
    if case.use_CPLEX == True:
        if (CPLEX()).available() == True:
            problem.solve(CPLEX())
        else:
            sys.stdout = old_stdout
            print('CPLEX is not available or not correctly installed')
            exit()
    else:
        problem.solve()
    if case.log_level < 2:
        sys.stdout = old_stdout

    end = time.time()
    print(Style.NORMAL, end='')
    
    if case.log_level >= 1:
        print('Computation time :\t\t{0:.3f} s'.format(end-start))    
        print("Status :\t\t\t", LpStatus[problem.status])
        print("Objective : \t\t\t{:.3f} TWh".format(value(problem.objective)/1e12))

    
    case.set_opti(e_exch, e_import, gas_exch, gas_import, cell_variables)
    case.compute_values()

    #return [e_exch, e_import, gas_exch, gas_import, problem, cell_variables,building]

    
def lp_problem_Pulp(cells, grid, policies, size, dt, use_CPLEX):
    return [grid_e_exch, grid_e_import, grid_gas_exch, grid_gas_import, lp, c_v, E_prod_tot, E_cons_tot, E_re_prod_tot, E_imp_tot, E_prod, E_re_prod, E_coal, E_nuc, E_ccgt, E_runriver, E_dam, E_pv, E_wt, ptg_use, E_imp, E_cons, curtailment]

