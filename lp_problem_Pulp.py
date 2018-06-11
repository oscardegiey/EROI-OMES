     
from pulp import *
import numpy as np
import time
    
    
def lp_problem_Pulp(cells, grid, policies, size, dt):
    
    start = time.time()
    lp = LpProblem('LP solver', LpMinimize)
    
    
    n = len(cells)
    
    re_share = policies[0]
    imp_share = policies[1]
    
    curtailment = {}
    
    E_invested = 0
    E_prod_tot = 0
    E_re_prod_tot = 0
    E_imp_tot = 0
    E_cons_tot = 0
    e_invested_battOff = 0
    e_invested_battOn = 0
    
    # Cell total prod
    E_prod = {}
    E_re_prod = {}
    E_coal = {}
    E_nuc = {}
    E_ccgt = {}
    E_runriver = {}
    E_dam = {}
    E_pv = {}
    E_wt = {}
    ptg_use = {}
        
    # Cell total import and consumption
    E_imp = {}
    E_cons = {}
    
    storage_weight = 0.5
    
    
    """ Grid variables initialisation """
    
    grid_e_exch = {}
    grid_gas_exch = {}
    grid_e_import = {}
    grid_gas_import = {}
    for cell1 in cells:
        grid_e_exch[cell1.name] = {}
        grid_gas_exch[cell1.name] = {}
        grid_e_import[cell1.name] = [0 for i in range(size)]
        grid_gas_import[cell1.name] = LpVariable('grid_gas_import' + cell1.name, 0, None)
        for cell2 in cells:
            grid_e_exch[cell1.name][cell2.name] = [0 for i in range(size)]
            grid_gas_exch[cell1.name][cell2.name]= [0 for i in range(size)]
    
    for t in range(size):
    
        # Electricity and gas exchange between cells (variable dictionnary)
        for cell1 in cells:
            for cell2 in cells:
            
                if cell1.name == cell2.name:
                    continue
            # No connection existing between the two cells
                elif grid.e_exch_max[cell1.name][cell2.name] == 0:
                    continue            
            # Cell1 variables dictionnary not yet created
                else:
                    grid_e_exch[cell1.name][cell2.name][t] = LpVariable('grid_e_exch_'+cell1.name+'_'+cell2.name+'_'+str(t), 0, grid.e_exch_max[cell1.name][cell2.name])            
                    
        for cell1 in cells:
            for cell2 in cells:
            
                if cell1.name == cell2.name:
                    continue
            # No connection existing between the two cells
                elif grid.gas_exch_max[cell1.name][cell2.name] == 0:
                    continue
            
            # Cell1 variables dictionnary not yet created
                else:
                    grid_gas_exch[cell1.name][cell2.name][t] = LpVariable('grid_gas_exch_'+cell1.name+'_'+cell2.name+'_'+str(t), -grid.gas_exch_max[cell1.name][cell2.name], grid.gas_exch_max[cell1.name][cell2.name])		
                 
        # Electricity and gas import from outside of the system in each cell
        for cell in cells:
            grid_e_import[cell.name][t] = LpVariable('grid_e_import' + cell.name+'_'+str(t), 0, grid.e_imp_max[cell.name])
    
    
    """ Cells variables initialisation """
    
    #Global variables dico
    c_v = {} 
    
    for cell in cells:
        #Cell variables dico
        v = {} 
        
        #Nuclear 
        #v['p_nuc'] =  LpVariable(cell.name+'p_nuc', 0, None)
        v['p_nuc'] =cell.nuclear.p_ut
        
        #Coal
        v['p_coal'] = [0 for i in range(size)]
        
        #CCGT
        v['p_ccgt'] = [0 for i in range(size)]
        
        #DAM
        v['p_out_dam'] = [0 for i in range(size)]
        v['p_in_dam'] = [0 for i in range(size)]
        v['e_dam'] = [0 for i in range(size)]
        
        #PV
        for park in cell.pv_parks:
            v[park.name+'_surf_pv'] = LpVariable(cell.name+park.name+'_surf_pv', park.minM2, park.maxM2) 
        #WT
        for park in cell.wt_parks: 
            v[park.name+'_n_wt'] = LpVariable(cell.name+park.name+'_n_wt', park.n_min, park.n_max)
        #RunRiver
        v['p_gen_runriver'] = [0 for i in range(size)]
        
        #PHES
        v['p_in_phes']= [0 for i in range(size)]
        v['p_out_phes']= [0 for i in range(size)]
        #v['p_pump_phes']= [0 for i in range(size)]
        v['e_phes']= [0 for i in range(size)]
        
        #Batteries
        for park in cell.bat_parks:
            v[park.type+'_capacity_bat'] = LpVariable(cell.name+park.type+'_capacity_bat'+'_'+str(t), park.Esizemin, park.Esizemax) 
            v[park.type+'_p_in_bat']= [0 for i in range(size)]
            v[park.type+'_p_out_bat']=[0 for i in range(size)]
            v[park.type+'_e_bat']= [0 for i in range(size)]
        
        #PtG
        v['p_in_ptg'] = [0 for i in range(size)]
        if cell.ptg.Pmax == 0:
            v['p_inst_ptg'] = 0	    
        else:
            v['p_inst_ptg'] = LpVariable(cell.name+'p_inst_ptg'+'_'+str(t), cell.ptg.Pmin, cell.ptg.Pmax) 
        
        #Gas
        if cell.gas.Esizemax == 0:
            v['capacity_gas'] = 0
        else:
            v['capacity_gas'] = LpVariable(cell.name+'capacity_gas'+'_'+str(t), cell.gas.Esizemin, cell.gas.Esizemax) 
        v['p_out_gas'] = [0 for i in range(size)]
        v['e_gas'] = [0 for i in range(size)]
        
        for t in range(size):
        #Producers
            
            #Coal
            if cell.coal.Pmax != 0:
                v['p_coal'][t] = LpVariable(cell.name+'p_coal'+'_'+str(t), cell.coal.Pmin, cell.coal.Pmax)
            
            #Ccgt
            if cell.ccgt.Pmax != 0:
                v['p_ccgt'][t] = LpVariable(cell.name+'p_ccgt'+'_'+str(t), cell.ccgt.Pmin, cell.ccgt.Pmax)        
                
            #dam
            if cell.dam.p_inst != 0:
                v['p_out_dam'][t] = LpVariable(cell.name+'p_out_dam'+'_'+str(t), 0, cell.dam.p_inst)
                v['p_in_dam'][t] = LpVariable(cell.name+'p_pin_dam'+'_'+str(t), -cell.dam.local_in[t], 0) #Negative
                v['e_dam'][t] = LpVariable(cell.name+'e_dam'+'_'+str(t), 0, cell.dam.capacity)
                
            #RunRiver
            if cell.runriver.p_inst != 0:
                v['p_gen_runriver'][t] = LpVariable(cell.name+'p_gen_runriver'+'_'+str(t), 0, cell.runriver.p_inst)
            
        #Storage
            
            #PHES
            if cell.phes.capacity != 0:
                v['p_in_phes'][t]= LpVariable(cell.name+'p_in_phes'+'_'+str(t), -cell.phes.p_inst_in, 0) #Negative
                v['p_out_phes'][t]=LpVariable(cell.name+'p_out_phes'+'_'+str(t), 0,cell.phes.p_inst_out)
                #v['p_pump_phes'][t]=LpVariable(cell.name+'p_pump_phes'+'_'+str(t), -cell.phes.p_inst_pump, 0) 
                v['e_phes'][t]= LpVariable(cell.name+'e_phes'+'_'+str(t), 0, cell.phes.capacity) 
            
                E_invested += -v['p_in_phes'][t]*dt/(cell.phes.ESOI)
            
            #Batteries 
            for park in cell.bat_parks:
                v[park.type+'_p_in_bat'][t]= LpVariable(cell.name+park.type+'_p_in_bat'+'_'+str(t),None, 0)   # Negative 
                v[park.type+'_p_out_bat'][t]=LpVariable(cell.name+park.type+'_p_out_bat'+'_'+str(t), 0, None)
                v[park.type+'_e_bat'][t]= LpVariable(cell.name+park.type+'_e_bat'+'_'+str(t), 0, None)
            
                E_invested += v[park.type+'_p_in_bat'][t]*park.eps_ut*(1-storage_weight)    
                if cell.name == 'Belgium_onshore':
                    e_invested_battOn += v[park.type+'_p_in_bat'][t]*park.eps_ut*(1-storage_weight)
                else:
                    e_invested_battOff += v[park.type+'_p_in_bat'][t]*park.eps_ut*(1-storage_weight)
                        
            #PtG
            if cell.ptg.Pmax != 0:

                v['p_in_ptg'][t] = LpVariable(cell.name+'p_in_ptg'+'_'+str(t),None , 0)  # Negative!!
                
                E_invested += v['p_in_ptg'][t]*cell.ptg.eps_ut*(1-storage_weight)    				
            
        #Gas
            if cell.gas.Esizemax != 0:
                v['e_gas'][t] = LpVariable(cell.name+'e_gas'+'_'+str(t), 0, None)
                
            
        # Store the variables in the variable dictionnary
        c_v[cell.name] = v
     
        
    """ Objective Function """
    
    for cell in cells:  
           
        for park in cell.pv_parks:
            E_invested += c_v[cell.name][park.name+'_surf_pv']*park.costpermsq*(size*dt)/park.lifetime
            
        for park in cell.wt_parks:
            E_invested += c_v[cell.name][park.name+'_n_wt']*park.costperWTinst*(size*dt)/park.lifetime
            
        for park in cell.bat_parks:
            E_invested += c_v[cell.name][park.type+'_capacity_bat']*park.eps_inst*(storage_weight)*(size*dt)/park.lifetime
            if cell.name == 'Belgium_onshore':
                e_invested_battOn += c_v[cell.name][park.type+'_capacity_bat']*park.eps_inst*(storage_weight)*(size*dt)/park.lifetime
            else:
                e_invested_battOff += c_v[cell.name][park.type+'_capacity_bat']*park.eps_inst*(storage_weight)*(size*dt)/park.lifetime
            
        E_invested += c_v[cell.name]['capacity_gas']*cell.gas.eps*(size*dt)/cell.gas.lifetime 
        
        E_invested += c_v[cell.name]['p_inst_ptg']*cell.ptg.eps_inst*(size*dt)/cell.ptg.lifetime*(storage_weight) 
        
    
    #objective
    lp += E_invested, 'E_invested'	
    
    
    """ Constraints """
    
    for cell in cells:
        # Cell total prod
        E_prod[cell.name] = 0
        E_re_prod[cell.name] = 0
        E_coal[cell.name] = 0
        E_nuc[cell.name] = 0
        E_ccgt[cell.name] = 0
        E_runriver[cell.name] = 0
        E_dam[cell.name] = 0
        E_pv[cell.name] = 0
        E_wt[cell.name] = 0
        ptg_use[cell.name] = 0
        
            
        # Cell total import and consumption
        E_imp[cell.name] = 0
        E_cons[cell.name] = 0
        
        curtailment[cell.name] = 0
        
        # PtG utilisation
        prod = np.zeros(size)
        for t in range(size):
            
            # Power balance 
            lp += (sum(c_v[cell.name][park.name+'_surf_pv']*park.prod_unit[t] for park in cell.pv_parks) + \
            sum(c_v[cell.name][park.name+'_n_wt']*park.prod_unit[t] for park in cell.wt_parks) + \
            c_v[cell.name]['p_ccgt'][t] + c_v[cell.name]['p_coal'][t] + c_v[cell.name]['p_nuc'] + \
            sum(c_v[cell.name][park.type+'_p_in_bat'][t] for park in cell.bat_parks) + \
            sum(c_v[cell.name][park.type+'_p_out_bat'][t] for park in cell.bat_parks) + \
            c_v[cell.name]['p_in_phes'][t] + c_v[cell.name]['p_out_phes'][t] + \
            c_v[cell.name]['p_out_dam'][t] + c_v[cell.name]['p_gen_runriver'][t] + \
            c_v[cell.name]['p_in_ptg'][t] - \
            sum(grid_e_exch[cell.name][neighbour.name][t] for neighbour in cells) + \
            sum(grid_e_exch[neighbour.name][cell.name][t] for neighbour in cells) + \
            grid_e_import[cell.name][t]  >= cell.cons[t])
            
            # Power exchange in both direction
            for cell2 in cells:
                lp += (grid_gas_exch[cell.name][cell2.name][t] == -grid_gas_exch[cell2.name][cell.name][t])

            # Storage balance

            # Batteries
            if t == 0:
                # Cyclic condition
                for park in cell.bat_parks:
                    lp += (c_v[cell.name][park.type + '_e_bat'][t] == c_v[cell.name][park.type + '_e_bat'][size-1]*pow((1-park.leak),dt) - \
                    c_v[cell.name][park.type + '_p_in_bat'][size-1]*park.efficiency_in*dt - \
                    c_v[cell.name][park.type + '_p_out_bat'][size-1]*1/park.efficiency_out*dt)
            else: 
                for park in cell.bat_parks:
                    lp += (c_v[cell.name][park.type + '_e_bat'][t] == c_v[cell.name][park.type + '_e_bat'][t-1]*pow((1-park.leak),dt) - \
                    c_v[cell.name][park.type + '_p_in_bat'][t-1]*park.efficiency_in*dt - \
                    c_v[cell.name][park.type + '_p_out_bat'][t-1]*1/park.efficiency_out*dt)

            # PHES
            if t == 0:
                # Cyclic condition
                lp += (c_v[cell.name]['e_phes'][t] == c_v[cell.name]['e_phes'][size-1] - \
                c_v[cell.name]['p_in_phes'][size-1]*cell.phes.efficiency_in*dt - \
                c_v[cell.name]['p_out_phes'][size-1]*1/cell.phes.efficiency_out*dt) 
            else:  
                lp += (c_v[cell.name]['e_phes'][t] == c_v[cell.name]['e_phes'][t-1] - \
                c_v[cell.name]['p_in_phes'][t-1]*cell.phes.efficiency_in*dt - \
                c_v[cell.name]['p_out_phes'][t-1]*1/cell.phes.efficiency_out*dt) 

            # Dam
            if t == 0:
                # Cyclic condition
                lp += (c_v[cell.name]['e_dam'][t] == c_v[cell.name]['e_dam'][size-1] - \
                c_v[cell.name]['p_in_dam'][size-1]*dt - \
                c_v[cell.name]['p_out_dam'][size-1]*1/cell.dam.efficiency_out*dt) 
            else: 
                lp += (c_v[cell.name]['e_dam'][t] == c_v[cell.name]['e_dam'][t-1] - \
                c_v[cell.name]['p_in_dam'][t-1]*dt - \
                c_v[cell.name]['p_out_dam'][t-1]*1/cell.dam.efficiency_out*dt) 
            
            
            # Gas
            if t == 0:
                # Cyclic condition
                lp += (c_v[cell.name]['e_gas'][t] == c_v[cell.name]['e_gas'][size-1] - \
                c_v[cell.name]['p_in_ptg'][size-1]*cell.ptg.efficiency_in*dt - \
                c_v[cell.name]['p_ccgt'][size-1]*1/cell.ccgt.efficiency*dt - \
                sum(grid_gas_exch[cell.name][neighbour.name][size-1]*dt for neighbour in cells) + \
                grid_gas_import[cell.name]*dt)
            else: 
                lp += (c_v[cell.name]['e_gas'][t] == c_v[cell.name]['e_gas'][t-1] - \
                c_v[cell.name]['p_in_ptg'][t-1]*cell.ptg.efficiency_in*dt - \
                c_v[cell.name]['p_ccgt'][t-1]*1/cell.ccgt.efficiency*dt - \
                sum(grid_gas_exch[cell.name][neighbour.name][t-1]*dt for neighbour in cells) + \
                grid_gas_import[cell.name]*dt)
                
            
            #Limits            
            
            #Batteries
            for park in cell.bat_parks:
                lp += (-c_v[cell.name][park.type+'_capacity_bat']*park.alpha_in<= c_v[cell.name][park.type+'_p_in_bat'][t]) # Power inlet
                lp += (c_v[cell.name][park.type+'_p_out_bat'][t]<=c_v[cell.name][park.type+'_capacity_bat']*park.alpha_out) # Power outlet
                lp += (c_v[cell.name][park.type+'_capacity_bat']*(1-park.DoD)<=c_v[cell.name][park.type+'_e_bat'][t]<=c_v[cell.name][park.type+'_capacity_bat']) # Energy stored
            
            
            # #DAM
            # if cell.dam.p_inst != 0:
            #     lp += (-cell.dam.local_in[t] <= c_v[cell.name]['p_in_dam'][t])
            
            #RunRiver
            if cell.runriver.p_inst != 0:
                lp += (c_v[cell.name]['p_gen_runriver'][t]<= cell.runriver.local_in[t]*cell.runriver.efficiency) #TODO: Est-ce que l'efficiency est vraiment utile?
            
            #PtG
            if cell.ptg.Pmax != 0:
                lp += (-c_v[cell.name]['p_inst_ptg'] <= c_v[cell.name]['p_in_ptg'][t])
                
            
            #GAS
            if cell.gas.Esizemax != 0:
                lp += (c_v[cell.name]['e_gas'][t]<=c_v[cell.name]['capacity_gas'])
                
                
            # Cell production
            # RE
            for park in cell.pv_parks:
                E_pv[cell.name] += c_v[cell.name][park.name+'_surf_pv']*park.prod_unit[t]*dt
            
            for park in cell.wt_parks:
                E_wt[cell.name] += c_v[cell.name][park.name+'_n_wt']*park.prod_unit[t]*dt
                
            #Dam
            E_dam[cell.name] += c_v[cell.name]['p_out_dam'][t]*dt
            
            #Runriver
            E_runriver[cell.name] += c_v[cell.name]['p_gen_runriver'][t]*dt
            
            # Non-RE
            #Ccgt
            E_ccgt[cell.name] += c_v[cell.name]['p_ccgt'][t]*dt
            
            
            #Nuc
            E_nuc[cell.name] += c_v[cell.name]['p_nuc']*dt
            
            #Coal
            E_coal[cell.name] += c_v[cell.name]['p_coal'][t]*dt
            
            # Cell import
            E_imp[cell.name]  += grid_e_import[cell.name][t]*dt
            
            # Cell consumption
            E_cons[cell.name]  += cell.cons[t]*dt
            
            # PtG utilisation
            ptg_use[cell.name] += c_v[cell.name]['p_in_ptg'][t]*dt
            
            
            
            curtailment[cell.name] += (sum(c_v[cell.name][park.name+'_surf_pv']*park.prod_unit[t] for park in cell.pv_parks) + \
            sum(c_v[cell.name][park.name+'_n_wt']*park.prod_unit[t] for park in cell.wt_parks) + \
            c_v[cell.name]['p_ccgt'][t] + c_v[cell.name]['p_coal'][t] + c_v[cell.name]['p_nuc'] + \
            sum(c_v[cell.name][park.type+'_p_in_bat'][t] for park in cell.bat_parks) + \
            sum(c_v[cell.name][park.type+'_p_out_bat'][t] for park in cell.bat_parks) + \
            c_v[cell.name]['p_in_phes'][t] + c_v[cell.name]['p_out_phes'][t] + \
            c_v[cell.name]['p_out_dam'][t] + c_v[cell.name]['p_gen_runriver'][t] + \
            c_v[cell.name]['p_in_ptg'][t] - \
            sum(grid_e_exch[cell.name][neighbour.name][t] for neighbour in cells) + \
            sum(grid_e_exch[neighbour.name][cell.name][t] for neighbour in cells) + \
            grid_e_import[cell.name][t]  - cell.cons[t])*dt
        
        
        
        E_re_prod[cell.name] = E_pv[cell.name]+E_wt[cell.name]+E_dam[cell.name]+E_runriver[cell.name] #- curtailment[cell.name]
        E_prod[cell.name] =  E_re_prod[cell.name]+grid_gas_import[cell.name]*cell.ccgt.efficiency*size*dt+E_nuc[cell.name]+E_coal[cell.name] 
        
        # Minimum RE share in cell production 
        lp += ((E_prod[cell.name]+E_imp[cell.name])*cell.re_share <= E_re_prod[cell.name])
        
        # Maximum cell import
        lp += (E_imp[cell.name] <= imp_share*E_cons[cell.name])
        
        # Nuc Share
        #lp+= (c_v[cell.name]['p_nuc']*size*dt<=(E_prod[cell.name]+E_imp[cell.name])*cell.nuc_share)
        
        # Total production
        E_prod_tot += E_prod[cell.name]
        E_re_prod_tot += E_re_prod[cell.name]
        
        #Total Import and consumption
        E_imp_tot += E_imp[cell.name]
        E_cons_tot += E_cons[cell.name]
        
        
        # PtG minimal use
        lp += (ptg_use[cell.name] <= -cell.ptg.alpha*c_v[cell.name]['p_inst_ptg']*size*dt) 
    
    
    # Minimum RE share in total production
    lp += ((E_prod_tot+E_imp_tot)*re_share <= E_re_prod_tot)
    
    # Maximum import
    lp += (E_imp_tot <= imp_share*E_cons_tot)
    
    end = time.time()
    print('Variables and constraints creation = ',end-start,'s')
             
                    
    
    lp.writeLP("France_Germany_Benelux.lp")
    
    #LpSolverDefault.msg = 1
    
    start = time.time()
    lp.solve(CPLEX())
    end = time.time()
    
    print('Computation time = ',end-start,'s')
    
    
    return [grid_e_exch, grid_e_import, grid_gas_exch, grid_gas_import, lp, c_v, E_prod_tot, E_cons_tot, E_re_prod_tot, E_imp_tot, E_prod, E_re_prod, E_coal, E_nuc, E_ccgt, E_runriver, E_dam, E_pv, E_wt, ptg_use, E_imp, E_cons, curtailment]

