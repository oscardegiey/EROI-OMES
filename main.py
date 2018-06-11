    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 09:42:23 2018

@author: thibaultmartinelle
"""


from get_France_Cell import *
from get_Germany_Benelux_Cell import *
from get_Italy_and_Alpine_states_Cell import *
from get_Iberian_Peninsula_Cell import *
from get_British_Isles_Cell import * 
from get_Scandinavia_Cell import *
from Class_Grid import *
from lp_problem_Pulp import *
from rebuild_problem import *
from display_data import *
import matplotlib.pyplot as plt
import numpy as np
from compute_val import *

plt.close('all')

from save_data import * 

plt.close('all')





dt = 128



size = int(8760/dt)
year = "2017"

""" Policies global system"""


P_RE = 1
P_imp = 0


policies = [P_RE, P_imp]


""" Cells initialisation """

re_share_fr = 0
imp_share_fr = 0.1
nuc_lf_fr = 0



re_share_de = 0
imp_share_de = 0.1
nuc_lf_de = 0


re_share_it = 0
imp_share_it = 0.1
nuc_lf_it = 0 

re_share_ib = 0
imp_share_ib = 0.1
nuc_lf_ib = 0 

re_share_bi = 0
imp_share_bi = 0.1
nuc_lf_bi = 0 

re_share_sc = 0
imp_share_sc = 0.1
nuc_lf_sc = 0 

France = get_France_Cell(re_share_fr, imp_share_fr, nuc_lf_fr, size, dt, year)
Germany_Benelux =  get_Germany_Benelux_Cell(re_share_de, imp_share_de, nuc_lf_de, size, dt, year)
Italy_and_Alpine_states = get_Italy_and_Alpine_states_Cell(re_share_it, imp_share_it, nuc_lf_it, size, dt, year)
Iberian_Peninsula = get_Iberian_Peninsula_Cell(re_share_ib,imp_share_ib,nuc_lf_ib,size,dt,year)
British_Isles = get_British_Isles_Cell(re_share_bi,imp_share_bi,nuc_lf_bi,size,dt,year)
Scandinavia = get_Scandinavia_Cell(re_share_sc,imp_share_sc,nuc_lf_sc,size,dt,year)

    

cells = [France, Germany_Benelux, Italy_and_Alpine_states, Iberian_Peninsula, British_Isles, Scandinavia]

# cells = [Germany_Benelux, France]



""" Import grid """

e_imp_max = {}
gas_imp_max = {}
e_imp_max['France'] = 0
e_imp_max['Germany_Benelux']= 0
e_imp_max['Italy_and_Alpine_states']= 0
e_imp_max['Iberian_Peninsula']= 0
e_imp_max['British_Isles']= 0
e_imp_max['Scandinavia']= 0
gas_imp_max['France']=0
gas_imp_max['Germany_Benelux']=0 
gas_imp_max['Italy_and_Alpine_states']= 0
gas_imp_max['Iberian_Peninsula']= 0
gas_imp_max['British_Isles']= 0
gas_imp_max['Scandinavia']= 0


# """ France Benelux case """
# e_imp_max['France'] = 2000e6 + 3100e6 + 2095e6
# e_imp_max['Germany_Benelux']= 2695e6 + 5500e6

 
""" Exchange grid """

e_exch_max = {
'France' : 
{'France' : 0, 'Germany_Benelux' : 6100e6*1, 'Italy_and_Alpine_states': 5775e6,                           'Iberian_Peninsula': 3000e6, 'British_Isles': 2000e6, 'Scandinavia': 0}, 
'Germany_Benelux' : 
{'France' : 5500e6*1, 'Germany_Benelux' : 0, 'Italy_and_Alpine_states': 3700e6, 'Iberian_Peninsula':0, 'British_Isles':0, 'Scandinavia': 2150e6}, 
'Italy_and_Alpine_states' : 
{'France' : 2095e6, 'Germany_Benelux' : 5500e6, 'Italy_and_Alpine_states': 0, 'Iberian_Peninsula': 0, 'British_Isles':0, 'Scandinavia': 0}, 
'Iberian_Peninsula': 
{'France' : 3100e6, 'Germany_Benelux' : 0, 'Italy_and_Alpine_states': 0, 'Iberian_Peninsula': 0, 'British_Isles': 0, 'Scandinavia': 0}, 
'British_Isles': 
{'France' : 2000e6, 'Germany_Benelux' : 0, 'Italy_and_Alpine_states': 0, 'Iberian_Peninsula': 0, 'British_Isles': 0, 'Scandinavia': 0},
'Scandinavia':
{'France' : 0, 'Germany_Benelux' : 2695e6, 'Italy_and_Alpine_states': 0, 'Iberian_Peninsula': 0, 'British_Isles': 0, 'Scandinavia': 0}}


gas_exch_max = {
'France' : 
{'Germany_Benelux' : 0, 'France' : 0, 'Italy_and_Alpine_states': 0, 'Iberian_Peninsula': 0, 'British_Isles': 0, 'Scandinavia': 0}, 
'Germany_Benelux' : 
{'France' : 0, 'Germany_Benelux' : 0, 'Italy_and_Alpine_states': 0, 'Iberian_Peninsula': 0,'British_Isles': 0, 'Scandinavia': 0}, 
'Italy_and_Alpine_states' : 
{'France' : 0, 'Germany_Benelux' : 0, 'Italy_and_Alpine_states': 0, 'Iberian_Peninsula': 0,'British_Isles': 0, 'Scandinavia': 0}, 
'Iberian_Peninsula' : 
{'France' : 0, 'Germany_Benelux' : 0, 'Italy_and_Alpine_states': 0, 'Iberian_Peninsula': 0, 'British_Isles': 0, 'Scandinavia': 0},
'British_Isles' :
{'France' : 0, 'Germany_Benelux' : 0, 'Italy_and_Alpine_states': 0, 'Iberian_Peninsula': 0, 'British_Isles': 0, 'Scandinavia': 0},
'Scandinavia':
{'France' : 0, 'Germany_Benelux' : 0, 'Italy_and_Alpine_states': 0, 'Iberian_Peninsula': 0, 'British_Isles': 0, 'Scandinavia': 0}}


grid = Grid(e_imp_max, e_exch_max, gas_imp_max, gas_exch_max)


"""Lp probelm """
print('Launch LP Problem')
[grid_e_exch, grid_e_import, grid_gas_exch, grid_gas_import, lp, c_v, E_prod_tot, E_cons_tot, E_re_prod_tot, E_imp_tot, E_prod, E_re_prod, E_coal, E_nuc, E_ccgt, E_runriver, E_dam, E_pv, E_wt, ptg_use, E_imp, E_cons, curtailment] = lp_problem_Pulp(cells, grid, policies, size, dt)

# The status of the solution is printed to the screen
print("Status:", LpStatus[lp.status])

# The optimised objective function value is printed to the screen
print("Objective = ", value(lp.objective))

""" Rebuilding the cells with the optimal values """
rebuild_problem(cells, grid, grid_e_exch, grid_e_import, grid_gas_exch, grid_gas_import,c_v,size,dt)


""" Saving the datas """
save_data(cells,grid,'eu_RE'+ str(P_RE) + '_imp'+ str(P_imp) +'_dt'+str(dt)+'.pkl')
#save_data(cells,grid,'FRGE1.1'+'.pkl')

# Print results in terminal

[EROI_PV, EROI_WT, Share,LF, E_prod, E_invested, E_invest_tot, E_invested_0] = compute_val(cells, grid)


print_data(cells,grid,Share,LF, E_prod,  E_invested, E_invest_tot, E_invested_0) #TODO: make automatic EROI print

# Plot data
ind = 0
#TODO: on reçoit un warning quand on ouvre full fenêtres
# ind = plot_cell_data(cells[0],size,dt,0)
# ind = plot_cell_data(cells[1],size,dt,ind)

for i in range(len(cells)):
    stack_plot_prod(cells,i,grid,size,dt)


plt.show()




