from pulp import *
import numpy as np

def rebuild_problem(cells, grid, grid_e_exch, grid_e_import, grid_gas_exch, grid_gas_import, variables_dico,size,dt):
    v = variables_dico
    
    i = 0
    
    gas_import = {}
    e_import = {}
    e_exch = {}
    gas_exch = {}
    
    for cell in cells:
        e_exch[cell.name] = {}
        gas_exch[cell.name] = {}
        
        """ Producers """
        
        for park in cell.pv_parks:
            surf_pv = value(v[cell.name][park.name+'_surf_pv'])
            p_gen = np.zeros(size)
            for t in range(size):
                p_gen[t] = surf_pv*park.prod_unit[t]
            park.set_opti(surf_pv,p_gen)           
        
        for park in cell.wt_parks:
            n_wt = value(v[cell.name][park.name+'_n_wt'])
            p_gen = np.zeros(size)
            for t in range(size):
                p_gen[t] = n_wt*park.prod_unit[t]
            park.set_opti(n_wt,p_gen)
        
            
        ccgt = np.zeros(size)
        coal = np.zeros(size)
        p_out_dam = np.zeros(size)
        e_dam = np.zeros(size)
        p_gen_runriver = np.zeros(size)
        for t in range(size):
            ccgt[t] = value(v[cell.name]['p_ccgt'][t])
            coal[t] = value(v[cell.name]['p_coal'][t])
            p_out_dam[t] = value(v[cell.name]['p_out_dam'][t])
            e_dam[t] = value(v[cell.name]['e_dam'][t])
            p_gen_runriver[t] = value(v[cell.name]['p_gen_runriver'][t])
        cell.ccgt.set_opti(ccgt)
        cell.coal.set_opti(coal)
        cell.dam.set_opti(p_out_dam,e_dam)
        cell.runriver.set_opti(p_gen_runriver)
        
        nuclear = value(v[cell.name]['p_nuc'])
        cell.nuclear.set_opti(nuclear)
        
        if nuclear>cell.nuclear.p_inst:
            print('Nuclear production exceed nuclear installed capacity. Change nuclear share.')
                
        
        
        """ Storages """
        
        #Batteries
        for park in cell.bat_parks:
            capacity = value(v[cell.name][park.type+'_capacity_bat'])
            energy = np.zeros(size)
            p_in = np.zeros(size)
            p_out = np.zeros(size)
            for t in range(size):
                energy[t] = value(v[cell.name][park.type+'_e_bat'][t])
                p_in[t] = value(v[cell.name][park.type+'_p_in_bat'][t])
                p_out[t] = value(v[cell.name][park.type+'_p_out_bat'][t])
            park.set_opti(capacity,energy,p_in,p_out)
        
        # PHES
        energy = np.zeros(size)
        p_in = np.zeros(size)
        p_out = np.zeros(size)
        #p_pump = np.zeros(size)
        for t in range(size):
            energy[t] = value(v[cell.name]['e_phes'][t])
            p_in[t] = value(v[cell.name]['p_in_phes'][t])
            p_out[t] = value(v[cell.name]['p_out_phes'][t])
           # p_pump[t] = value(v[cell.name]['p_pump_phes'][t])
        cell.phes.set_opti(energy,p_in,p_out)
        
        # PtG
        p_inst = value(v[cell.name]['p_inst_ptg'])
        p_in = np.zeros(size)
        for t in range(size):
            p_in[t] = value(v[cell.name]['p_in_ptg'][t])
        cell.ptg.set_opti(p_inst,p_in)
        
        #Gas
        capacity = value(v[cell.name]['capacity_gas'])
        energy = np.zeros(size)
        for t in range(size):
            energy[t] = value(v[cell.name]['e_gas'][t])
        cell.gas.set_opti(capacity,energy)
        
        
        """ Grid """
        
        # Import
        gas_import[cell.name]= value(grid_gas_import[cell.name])
        e_import[cell.name] = np.zeros(size)
        for t in range(size):
            e_import[cell.name][t] = value(grid_e_import[cell.name][t])
            
        # Exchange
        
        for neighbour in cells:
            e_exch[cell.name][neighbour.name] = np.zeros(size)
            gas_exch[cell.name][neighbour.name] = np.zeros(size)
            for t in range(size):
                e_exch[cell.name][neighbour.name][t] = value(grid_e_exch[cell.name][neighbour.name][t])
                gas_exch[cell.name][neighbour.name][t] = value(grid_gas_exch[cell.name][neighbour.name][t])
            
    grid.set_opti(e_import, e_exch, gas_import, gas_exch)
            
            

        