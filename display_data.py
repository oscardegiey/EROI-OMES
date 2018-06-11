

from matplotlib.pyplot import *
from pulp import *
import numpy as np
from prettytable import PrettyTable # /!\ Need to install prettytable package with pip
import sys

def plot_cell_data(cell,size,dt,ind):
    time = np.linspace(1,size*dt,size)
    
    i = ind
    
    """ Conso """
    figure(i)
    plot(time,cell.cons)
    title('Electricity consumption '+cell.name)
    xlabel('Time [h]')
    ylabel('Power [W]')
    i +=1
    
    
    """ Producers """
    
    for park in cell.pv_parks:
        figure(i)
        plot(time,park.p_gen)
        title('PV production in cell '+cell.name+', park '+park.name)
        xlabel('Time [h]')
        ylabel('Power [W]')
        i +=1
        
    for park in cell.wt_parks:
        figure(i)
        plot(time,park.p_gen)
        title('Wind production in cell '+cell.name+', park '+park.name)
        xlabel('Time [h]')
        ylabel('Power [W]')
        i +=1
    
        figure(i)
        plot(time,cell.ccgt.p_gen)
        title('CCGT production in cell '+cell.name)
        xlabel('Time [h]')
        ylabel('Power [W]')
        i +=1
        
        figure(i)
        plot(time,cell.dam.p_out)
        title('Dam production in cell '+cell.name)
        xlabel('Time [h]')
        ylabel('Power [W]')
        i +=1
        
        figure(i)
        plot(time,cell.dam.energy)
        title('Dam storage in cell '+cell.name)
        xlabel('Time [h]')
        ylabel('Energy [Wh]')
        i +=1
        
        figure(i)
        plot(time,cell.runriver.p_gen)
        title('RunRiver production in cell '+cell.name)
        xlabel('Time [h]')
        ylabel('Power [W]')
        i +=1
        
        figure(i)
        plot(time,cell.coal.p_gen)
        title('Coal production in cell '+cell.name)
        xlabel('Time [h]')
        ylabel('Power [W]')
        i +=1
        
    """ Storages """ 
    # Batteris
    for park in cell.bat_parks:
        figure(i)
        plot(time,park.energy)
        title('Battery storage in cell '+cell.name+', park '+park.type)
        xlabel('Time [h]')
        ylabel('Energy [Wh]')
        
        figure(i+1)
        plot(time,park.p_in)
        title('Battery input power in cell '+cell.name+', park '+park.type)
        xlabel('Time [h]')
        ylabel('Power [W]')
        
        figure(i+2)
        plot(time,park.p_out)
        title('Battery output power in cell '+cell.name+', park '+park.type)
        xlabel('Time [h]')
        ylabel('Power [W]')
        i +=3
    # PHES  
    figure(i)
    plot(time,cell.phes.energy)
    title('PHES storage in cell '+cell.name)
    xlabel('Time [h]')
    ylabel('Energy [Wh]')
    
    figure(i+1)
    plot(time,cell.phes.p_in)
    title('PHES input power in cell '+cell.name + '(local + p_pump - spillage)')
    xlabel('Time [h]')
    ylabel('Power [W]')
    
    figure(i+2)
    plot(time,cell.phes.p_out)
    title('PHES output power in cell '+cell.name)
    xlabel('Time [h]')
    ylabel('Power [W]')
    
    # figure(i+3)
    # plot(time,cell.phes.p_pump)
    # title('PHES pumped power in cell '+cell.name)
    # xlabel('Time [h]')
    # ylabel('Power [W]')
    
    i +=3
    
    # PtG
    
    figure(i)
    plot(time,cell.ptg.p_in)
    title('PtG input power in cell '+cell.name)
    xlabel('Time [h]')
    ylabel('Power [W]')
    
    i += 1
    
    # Gas
    
    figure(i)
    plot(time,cell.gas.energy)
    title('Gas storage in cell '+cell.name)
    xlabel('Time [h]')
    ylabel('Energy [Wh]')
    
    i += 2
    
    
        
    return i
    
    
def stack_plot_prod(cells,cell_nbr,grid,size,dt):
    
    name_dic = {'Germany_Benelux': 'Germany&Benelux', 'France': 'France', 'Iberian_Peninsula': 'Iberian Peninsula', 'Italy_and_Alpine_states': 'Italy and Alpine states', 'British_Isles': 'British Isles', 'Scandinavia': 'Scandinavia'}
    
    
    cell = cells[cell_nbr]
    
    p_unit = 1e9
    
    time = np.linspace(1,size*dt,size)
    legend_prod = []
    color_prod = []
    
    prod = []
    
    # prod.append(cell.nuclear.p_gen)
    # legend_prod += ['Nuclear']
    # color_prod += ['mediumpurple']
    
    prod.append(cell.runriver.p_gen)
    legend_prod += ['Runriver']
    color_prod += ['#42FFFC']
    
    prod.append(cell.dam.p_out)
    legend_prod += ['Dam']
    color_prod += ['#0003c7']
    
    for park in cell.pv_parks:
        prod.append(park.p_gen)
        legend_prod.append('PV')
    color_prod += ['#ffe402']
    
    for park in cell.wt_parks:
        prod.append(park.p_gen)
    color_prod += ['#a3e2a1', '#00B54E']
    legend_prod += ['Wind onshore', 'Wind offshore']
    
    prod.append(cell.ccgt.p_gen)
    legend_prod += ['CCGT']
    color_prod += ['#ff7373']
     
    for park in cell.bat_parks:
        prod.append(park.p_out)
        legend_prod.append('Battery out')
    color_prod += ['#606060'] 
        
    prod.append(cell.phes.p_out)
    legend_prod.append('PHES p_out')
    color_prod += ['#8ce2ff']
    
    # prod.append(cell.coal.p_gen)
    # legend_prod += ['Coal']
    # color_prod += ['#7A3B00']
    
    
    e_exch_in = np.zeros(size)
    for neighbour in cells:
        e_exch_in = e_exch_in + grid.e_exch[neighbour.name][cell.name]
    
    prod.append(e_exch_in)
    legend_prod.append('Import (other cells)')
    color_prod += ['#FFBB87']
    
    # prod.append(grid.e_import[cell.name])
    # legend_prod += ['Import (out system)']
    # color_prod += ['#FFBB87'] #TODO change color
    
    
    
    
    """Loads and plots"""
    
    figure()
    polys = stackplot(time,np.array(prod)/p_unit,colors = color_prod)
    
    load = 0
    load += cell.cons/p_unit
    cons, = plot(time,load,'r-', linestyle = '-', zorder = 5)
    legend_prod.append('Consumption')
    
    load += -cell.ptg.p_in/p_unit
    ptg, = plot(time,load,color = '#ff7373',linestyle = '-', zorder = 4)
    legend_prod.append('PtG')
    
    load += -cell.bat_parks[0].p_in/p_unit
    bat_in, = plot(time,load,color = '#606060', linestyle = '-', zorder = 3)
    legend_prod.append('Battery in')
    
    load += -cell.phes.p_in/p_unit
    phes_in, = plot(time,load,'#2791f4', linestyle = '-', zorder = 2)
    legend_prod.append('PHES in')
    
    # e_exch_out = np.zeros(size)
    # for neighbour in cells:
    #     e_exch_out += grid.e_exch[cell.name][neighbour.name]
    # load += e_exch_out/p_unit
    # export, = plot(time,load,'k', linestyle = '-', zorder = 1)
    # legend_prod.append('Export')


    rectangles = []
    for poly in polys:
        rectangles.append(Rectangle((0,0),1,1, fc=poly.get_facecolor()[0]))
    rectangles += [cons, ptg, bat_in, phes_in]
    
    legend(rectangles,legend_prod, loc=9, bbox_to_anchor=(0.5, -0.17), ncol=3, fontsize = 12)
    
    title('Production in\n' + name_dic[cell.name])
    xlabel('Time [h]')
    ylabel('Production [GW]')
    
    subplots_adjust(bottom = 0.37, left = 0.06, right = 0.98)
    



def stack_plot_storage(cells,cell_nbr,grid,size,dt,ind):
    
    cell = cells[cell_nbr]
    storage = []
    legend_storage = []
    
    time = np.linspace(1,size*dt,size)
    
    e_exch_out = np.zeros(size)
    e_exch_in = np.zeros(size)
    # for neighbour in cells:
    #     e_exch_out = e_exch_out - grid.e_exch[neighbour.name][cell.name]
    #     e_exch_out = e_exch_out + grid.e_exch[cell.name][neighbour.name]
    #     e_exch_out[e_exch_out<0] = 0
    #     e_exch_in = e_exch_in + grid.e_exch[neighbour.name][cell.name]
    #     e_exch_in = e_exch_in - grid.e_exch[cell.name][neighbour.name]
    #     e_exch_in[e_exch_in<0] = 0
    
    for neighbour in cells:
        e_exch_out = e_exch_out + grid.e_exch[cell.name][neighbour.name]
        e_exch_in = e_exch_in + grid.e_exch[neighbour.name][cell.name]
    
    storage.append(e_exch_out)
    legend_storage.append('Electricity export to other cells')
    
    for park in cell.bat_parks:
        storage.append(-park.p_in)
        legend_storage.append('Battery ' + park.type + ' p_in')
    
    storage.append(-cell.phes.p_in)
    legend_storage.append('PHES p_in')
    storage.append(-cell.ptg.p_in)
    legend_storage.append('PtG p_in')
    
    excess = sum(park.surf_pv*park.prod_unit for park in cell.pv_parks) + \
            sum(park.n_wt*park.prod_unit for park in cell.wt_parks) + \
            cell.ccgt.p_gen + cell.coal.p_gen + cell.nuclear.p_gen + \
            sum(park.p_out for park in cell.bat_parks) + \
            cell.phes.p_out + \
            cell.dam.p_out + cell.runriver.p_gen + \
            e_exch_in + \
            grid.e_import[cell.name] - cell.cons
    
    """Storage"""
    
    figure(ind)
    polys = stackplot(time,storage)
    
    excess_prod, = plot(time,excess,'b-')
    legend_storage.append('Excess production')
    
    rectangles = []
    for poly in polys:
        rectangles.append(Rectangle((0,0),1,1, fc=poly.get_facecolor()[0]))
    rectangles.append(excess_prod)
    
    legend(rectangles,legend_storage)
    title('Storage input in cell ' + cell.name)
    
    return ind+1



def plot_grid_data(cells, grid, size, dt,ind):
    time = np.linspace(1,size*dt,size)
    
    """ Grid """
    i = ind
    
    for cell in cells:
    
        # Import 
        figure(i)
        plot(time,grid.e_import[cell.name])
        title('Electricity import in cell' + cell.name)
        xlabel('Time [h]')
        ylabel('Power [W]')
        
        i += 1
        #TODO: add condition to avoid to print graph for zero transfer (possible with max)
        
            
        
        for neighbour in cells:
            if cell.name == neighbour.name:
                continue
            
            figure(i)
            plot(time,grid.e_exch[cell.name][neighbour.name])
            title('Electricity exchange cell ' + cell.name + ' --> ' + neighbour.name) 
            xlabel('Time [h]')
            ylabel('Power [W]')
            
            figure(i+1)
            plot(time,grid.gas_exch[cell.name][neighbour.name])
            title('Gas exchange cell ' + cell.name + ' --> ' + neighbour.name) 
            xlabel('Time [h]')
            ylabel('Power [W]')
            
            i += 2
    
    return i
    
    

def print_data(cells,grid,Share,LF, E_prod,  E_invested, E_invest_tot, E_invested_0):
    E_prod_tot = 0
    E_cons_tot = 0
    E_RE = 0
    E_imp_tot = 0
    EROI_PV = {}
    EROI_WT = {}
    E_cons = {}
    for cell in cells:
        E_cons[cell.name] = sum(cell.cons)*cell.dt
        EROI_PV[cell.name]= {}
        EROI_WT[cell.name] = {}
        E_prod_tot += (sum(cell.nuclear.p_gen)+grid.gas_import[cell.name]*cell.ccgt.efficiency*cell.size+sum(cell.coal.p_gen)+sum(cell.dam.p_out)+sum(cell.runriver.p_gen))*cell.dt
        E_cons_tot += sum(cell.cons)*cell.dt
        E_RE += (sum(cell.dam.p_out)+sum(cell.runriver.p_gen))*cell.dt
        for t in range(cell.size):
            E_imp_tot += grid.e_import[cell.name][t]*cell.dt
        for park in cell.pv_parks:
            EROI_PV[cell.name][park.name]= sum(park.prod_unit*cell.dt)/park.costpermsq*park.lifetime/(cell.size*cell.dt)
            E_prod_tot += sum(park.p_gen)*cell.dt
            E_RE += sum(park.p_gen)*cell.dt
        for park in cell.wt_parks:
            EROI_WT[cell.name][park.name] = sum(park.prod_unit)*cell.dt/park.costperWTinst*park.lifetime/(cell.size*cell.dt)
            E_prod_tot += sum(park.p_gen)*cell.dt
            E_RE += sum(park.p_gen)*cell.dt
 
    
    t_pv = PrettyTable(['PV Park','Surface [m^2]','Bounds [m^2]'])
    t_wt = PrettyTable(['WT Park','# of WT','Bounds'])
    t_bat = PrettyTable(['Battery Park', 'Installed capacity [Wh]','Bounds'])
    t_gas_1 = PrettyTable(['Cell', 'Installed PtG capacity [W]','Installed Gas capacity [Wh]'])
    t_gas_2 = PrettyTable(['Cell', 'CCGT inctalled capacity [W]','Gas import [W]'])
    t_nuc = PrettyTable(['Cell', 'Nuclear installed capacity [W]','Intsalled nuclear 2018[W]'])
    t_e_invest = PrettyTable(['Cell', 'Energy investement in RE assets','Energy investement in storage', 'EROI' ])
    t_e_invest_new= PrettyTable(['Cell', 'Energy investement in RE assets-Already installed','Energy investement in storage'])
    
    t_EROI_PV= PrettyTable(['Park', 'EROI PV'])
    t_EROI_WT= PrettyTable(['Park', 'EROI WT'])
    
    
    t_LF = PrettyTable(['Cell', 'Dam LF', 'Runriver LF','PtG LF','Batteries LF', 'CCGT LF'])
    
    t_tot = PrettyTable([])
    
    t_tot.add_row(['Total energy production [Wh]', value(E_prod_tot)])
    t_tot.add_row(['Total energy consumption [Wh]', value(E_cons_tot)])
    t_tot.add_row(['Total RE production [Wh]', value(E_RE)])
    t_tot.add_row(['Total electicity import [Wh]', '%e'%value(E_imp_tot)])
    
    t_share = {}
    E_invested_0_tot = 0
    E_cons_tot = 0
    
    for cell in cells:
        t_share[cell.name] =  PrettyTable([cell.name,'share'])
       
        for share in Share[cell.name]:
            t_share[cell.name].add_row([share, Share[cell.name][share]])
            
        for park in cell.pv_parks:
                t_pv.add_row([park.name, '%e'%park.surf_pv,'['+ str('%e'%park.minM2) + ' : ' + str('%e'%park.maxM2) + ']'])
                
                
        for park in cell.wt_parks:
                t_wt.add_row([park.name, park.n_wt,'['+ str(park.n_min) + ' : ' + str(park.n_max) + ']'])
            
            # Batteries
        for park in cell.bat_parks:
                t_bat.add_row([cell.name + ' ' + park.type, '%e'%park.capacity, '['+ str('%e'%park.Esizemin) + ' : ' + str('%e'%park.Esizemax) + ']'])
                
            # PtG and Gas 
            
        t_gas_1.add_row([cell.name, '%e'%cell.ptg.p_inst, '%e'%cell.gas.capacity])
            
            # Gas import and CCGT
        t_gas_2.add_row([cell.name, '%e'%cell.ccgt.Pmax, '%e'%grid.gas_import[cell.name]])
            
            # Nuclear
        t_nuc.add_row([cell.name,'%e'%cell.nuclear.p_gen[0], '%e'%cell.nuclear.p_inst ])
            
            #EROI
        for park in cell.pv_parks:
                t_EROI_PV.add_row([park.name, EROI_PV[cell.name][park.name]])
                
        for park in cell.wt_parks:
                t_EROI_WT.add_row([park.name, EROI_WT[cell.name][park.name]])
            
        # Load factor
        t_LF.add_row([cell.name,'%e'%LF[cell.name]['Dam'], '%e'%LF[cell.name]['Runriver'], '%e'%LF[cell.name]['PtG'], '%e'%LF[cell.name]['Batt'], '%e'%LF[cell.name]['CCGT']])
        
        # Enegy invested
        t_e_invest.add_row([cell.name, E_invested['REassets'][cell.name], E_invested['Storage'][cell.name], sum(cell.cons)*cell.dt/(E_invested['REassets'][cell.name] + E_invested['Storage'][cell.name])])
        t_e_invest_new.add_row([cell.name, E_invested['REassets'][cell.name]-E_invested_0[cell.name], E_invested['Storage'][cell.name]])
        E_invested_0_tot += E_invested_0[cell.name]
        E_cons_tot += sum(cell.cons)*cell.dt
        
    # Production shares 
    
    t_tot2 =  PrettyTable([' ','share', 'Prod RE/cons tot'])
    t_tot2.add_row(['RE share tot', Share['tot']['RE'], E_RE/E_cons_tot])
    t_prod = PrettyTable(['Cell', 'Production', 'Production/Production_tot', 'Production/conso','Production RE/conso'])
    t_e_invest.add_row(['GLOBAL', E_invest_tot['REassets'], E_invest_tot['Storage'],E_cons_tot/(E_invest_tot['REassets'] + E_invest_tot['Storage']) ])
    t_e_invest_new.add_row(['GLOBAL', E_invest_tot['REassets']-E_invested_0_tot, E_invest_tot['Storage']])
    for cell in cells: 
        t_prod.add_row([cell.name, E_prod[cell.name], E_prod[cell.name]/E_prod['tot'],E_prod[cell.name]/E_cons[cell.name], Share[cell.name]['RE/cons']])
    print(t_pv)
    print(t_wt)
    print(t_bat)
    print(t_gas_1)
    print(t_gas_2)
    print(t_nuc)
    print(t_tot)
    print(t_EROI_PV)
    print(t_EROI_WT)
    print(t_LF)
    print(t_tot2)
    print(t_prod)
    print(t_e_invest)
    print(t_e_invest_new)
    
    for cell in cells:
        print(t_share[cell.name])
    
    
def plot_month_storage(cells,grid,size,dt,ind):

    
    """Data processing"""
    
    
    fig_n = ind
    month = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,0])*24/dt
    m = np.zeros(13)
    for i in range(12):
        m[i] = m[i-1]+month[i]
    m = np.int16(m)
    
    # Units
    e_ratio = 1e12
    
    for cell in cells:
        p_m = {}
        p_m['batt'] = {}
        p_m['ptg'] = {}
        p_m['phes'] = {}
        p_m['batt_tot'] = np.zeros(12)
        p_m['ptg_tot'] = np.zeros(12)
        p_m['phes_tot'] = np.zeros(12)
        p_m['cons'] = np.zeros(12)
        p_m['export'] = np.zeros(12)
        
        prod_assets = {'pv','wt_on','wt_off','import','hydro'}
        for i in range(12):
            prod_tot = 0
            
            # Production
            prod = {}
            prod['pv'] = cell.pv_parks[0].prod_unit[m[i-1]:m[i]]*cell.pv_parks[0].surf_pv
            prod['wt_on'] = cell.wt_parks[0].prod_unit[m[i-1]:m[i]]*cell.wt_parks[0].n_wt
            prod['wt_off'] = cell.wt_parks[1].prod_unit[m[i-1]:m[i]]*cell.wt_parks[1].n_wt
            prod['hydro'] = cell.dam.p_out[m[i-1]:m[i]] + \
            cell.runriver.p_gen[m[i-1]:m[i]]
            
            # Storage
            storage = {}
            storage['batt'] = cell.bat_parks[0].p_in[m[i-1]:m[i]]
            storage['ptg'] = cell.ptg.p_in[m[i-1]:m[i]]
            storage['phes'] = cell.phes.p_in[m[i-1]:m[i]]
            
            # Import
            imp_global = 0
            for neighbour in cells:
                imp_global += grid.e_exch[neighbour.name][cell.name][m[i-1]:m[i]] - \
                grid.e_exch[cell.name][neighbour.name][m[i-1]:m[i]]
                
            prod['import'] = imp_global.clip(0)
            p_m['export'][i] = sum(imp_global.clip(max=0))*dt/e_ratio
            
            prod_tot = (prod['pv'] + prod['wt_on'] + prod['wt_off'] + prod['hydro'] + prod['import'])
            
            # Consumption
            p_m['cons'][i] = sum(cell.cons[m[i-1]:m[i]])*dt/e_ratio
            

            p_m['batt'][i] = {}
            p_m['ptg'][i] = {}
            p_m['phes'][i] = {}
            p_m['batt_tot'][i] = 0
            p_m['ptg_tot'][i] = 0
            p_m['phes_tot'][i] = 0
            for producer in prod_assets:
                p_m['batt'][i][producer] = sum(prod[producer]/prod_tot*storage['batt'])*dt/e_ratio
                p_m['ptg'][i][producer] = sum(prod[producer]/prod_tot*storage['ptg'])*dt/e_ratio
                p_m['phes'][i][producer] = sum(prod[producer]/prod_tot*storage['batt'])*dt/e_ratio
                p_m['batt_tot'][i] += p_m['batt'][i][producer]
                p_m['ptg_tot'][i] += p_m['ptg'][i][producer]
                p_m['phes_tot'][i] += p_m['phes'][i][producer]
         
        
        """ Plots """
        figure(fig_n)
        
        
        N = 12
        ind = np.arange(N)    # the x locations for the groups
        width = 0.35       # the width of the bars: can also be len(x) sequence
        
        color_assets = {'pv':'#ffe402','wt_on':'#72c245', 'wt_off': '#00B54E', 
        'import': '#000000' ,'hydro': '#42FFFC' }
        
        # Storage 
        btm = 0
        p0 = bar(ind-width/2, p_m['batt_tot'], width, bottom = btm, color = '#A0A0A0', label = 'Batteries') #grey
        btm += p_m['batt_tot']
        p1 = bar(ind-width/2, p_m['ptg_tot'], width, bottom = btm, color = '#cc3300', label = 'PtG' ) #grey
        btm += p_m['ptg_tot']
        p2 = bar(ind-width/2, p_m['phes_tot'], width, bottom = btm, color = '#66B2FF', label = 'PHES') #grey
        btm += p_m['phes_tot']
        
        # Consumption and export
        btm = 0
        p3 = bar(ind, p_m['cons'], width, bottom = btm, color = '#8666c1', label = 'Consumption') #grey
        btm += p_m['cons']
        p3 = bar(ind, -p_m['export'], width, bottom = btm, color = '#000000', label = 'Export') #grey
        
        
        btm = 0
        prod = np.zeros(12)
        for producer in prod_assets:
            for i in range(12):
                prod[i] = p_m['batt'][i][producer]
            bar(ind+width/2, prod, width, bottom = btm, color = color_assets[producer], label = producer) #grey
            btm += prod
        
        for producer in prod_assets:
            for i in range(12):
                prod[i] = p_m['ptg'][i][producer]
            bar(ind+width/2, prod, width, bottom = btm, color = color_assets[producer]) #grey
            btm += prod
            
        for producer in prod_assets:
            for i in range(12):
                prod[i] = p_m['phes'][i][producer]
            bar(ind+width/2, prod, width, bottom = btm, color = color_assets[producer]) #grey
            btm += prod
        
        
        ylabel('Production [TWh]')
        title('Monthly energy stored in '+ cell.name)
        xticks(ind, ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
        #yticks(np.arange(0, 81, 10))
        
        legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=4)
        
        #legend(['Consumption','Battery in', 'PtG', 'PHES in', 'e_exch_out', 'PV', 'WT onshore', 'WT offshore', 'Dam', 'Run-of-river','PHES out', 'Nuclear', 'Coal', 'CCGT', 'Battery out', 'e_exch_in', 'e_import'],loc=9, bbox_to_anchor=(0.5, -0.1), ncol=4)
        
        subplots_adjust(bottom=0.25)
        
        savefig('img/'+ cell.name+'_storage_monthly.eps', format='eps', dpi=1000)
        
        fig_n += 1
                
        
def plot_month_gas(cells,grid,size,dt,ind):
    #rcParams.update({'font.size': 15})
    
    name = {'Germany_Benelux': 'Germany&Benelux', 'France': 'France', 'Iberian_Peninsula': 'Iberian Peninsula', 'Italy_and_Alpine_states': 'Italy and Alpine states', 'British_Isles': 'British Isles', 'Scandinavia': 'Scandinavia', 'phes': 'PHES', 'ptg': 'PtG', 'bat': 'Batteries'}
    
    """Data processing"""
    fig_n = ind
    month = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,0])*24/dt
    m = np.zeros(13)
    for i in range(12):
        m[i] = m[i-1]+month[i]
    m = np.int16(m)
    
    # Units
    e_ratio = 1e12
    
    p_m_cells = {}
    for cell in cells:
        
        p_m = {}
        p_m['wt_on'] = np.zeros(12)
        p_m['wt_off'] = np.zeros(12)
        p_m['pv'] = np.zeros(12)
        p_m['dam'] = np.zeros(12)
        p_m['runriver'] = np.zeros(12)
        p_m['phes_out'] = np.zeros(12)
        p_m['phes_in'] = np.zeros(12)
        p_m['nuc'] = np.zeros(12)
        p_m['coal'] = np.zeros(12)
        p_m['ccgt'] = np.zeros(12)
        p_m['bat_out'] = np.zeros(12)
        p_m['bat_in'] = np.zeros(12)
        p_m['ptg'] = np.zeros(12)
        p_m['e_exch_out'] = np.zeros(12)
        p_m['e_exch_in'] = np.zeros(12)
        p_m['e_imp'] = np.zeros(12)
        p_m['e_bat'] = np.zeros(12)
        p_m['e_ptg'] = np.zeros(12)
        p_m['e_phes'] = np.zeros(12)
        
        for i in range(12):
            
            
            # Conventional
            p_m['ccgt'][i] = sum(cell.ccgt.p_gen[m[i-1]:m[i]])/cell.ccgt.efficiency*dt/e_ratio
            
           
            # PtG
            p_m['ptg'][i] = sum(cell.ptg.p_in[m[i-1]:m[i]])*cell.ptg.efficiency_in*dt/e_ratio
            

            # Storage
            #p_m['e_ptg'][i] = np.mean(cell.gas.energy[m[i-1]:m[i]])/e_ratio
            p_m['e_ptg'][i] = cell.gas.energy[m[i]-1]/e_ratio
                    
        p_m_cells[cell.name] = p_m
        
        
        
        """ Plots """
        figure(fig_n)
        
        
        N = 12
        ind = np.arange(N)    # the x locations for the groups
        width = 0.35       # the width of the bars: can also be len(x) sequence
        
        # Gas storage
        btm = 0
        p1 = bar(ind, -p_m['ptg']-p_m['ccgt'], width, bottom = btm, color = '#cc3300', label = 'PtG in') #red
        btm += -p_m['ptg']-p_m['ccgt']
        


        # Energy stored
        plot(ind, p_m['e_ptg'], '--o', color='b', label = 'Gas storage')

        
        ylabel('Energy In/Out [TWh]')
        title('Monthly energy stored in Gas storage in '+ name[cell.name])
        xticks(ind, ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
        #yticks(np.arange(0, 81, 10))
        
        legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=4)
        
        subplots_adjust(bottom=0.2)
        
        savefig('img/'+ cell.name+'_storage_monthly.eps', format='eps', dpi=1000)
        
        fig_n += 1

def plot_month_bar(cells,grid,size,dt,ind):
    
    #rcParams.update({'font.size': 15})
    
    name = {'Germany_Benelux': 'Germany&Benelux', 'France': 'France', 'Iberian_Peninsula': 'Iberian Peninsula', 'Italy_and_Alpine_states': 'Italy and Alpine states', 'British_Isles': 'British Isles', 'Scandinavia': 'Scandinavia', 'phes': 'PHES', 'ptg': 'PtG', 'bat': 'Batteries'}
    
    """Data processing"""
    fig_n = ind
    month = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,0])*24/dt
    m = np.zeros(13)
    for i in range(12):
        m[i] = m[i-1]+month[i]
    m = np.int16(m)
    
    # Units
    e_ratio = 1e12
    
    p_m_cells = {}
    for cell in cells:
        
        p_m = {}
        p_m['wt_on'] = np.zeros(12)
        p_m['wt_off'] = np.zeros(12)
        p_m['pv'] = np.zeros(12)
        p_m['dam'] = np.zeros(12)
        p_m['runriver'] = np.zeros(12)
        p_m['phes_out'] = np.zeros(12)
        p_m['phes_in'] = np.zeros(12)
        p_m['nuc'] = np.zeros(12)
        p_m['coal'] = np.zeros(12)
        p_m['ccgt'] = np.zeros(12)
        p_m['bat_out'] = np.zeros(12)
        p_m['bat_in'] = np.zeros(12)
        p_m['ptg'] = np.zeros(12)
        p_m['e_exch_out'] = np.zeros(12)
        p_m['e_exch_in'] = np.zeros(12)
        p_m['e_imp'] = np.zeros(12)
        p_m['cons'] = np.zeros(12)
        
        for i in range(12):
            
                       
            # WT
            p_m['wt_on'][i] = sum(cell.wt_parks[0].prod_unit[m[i-1]:m[i]] * cell.wt_parks[0].n_wt)*dt/e_ratio
            p_m['wt_off'][i] = sum(cell.wt_parks[1].prod_unit[m[i-1]:m[i]] * cell.wt_parks[1].n_wt)*dt/e_ratio
            
            # PV
            p_m['pv'][i] = sum(cell.pv_parks[0].prod_unit[m[i-1]:m[i]] * cell.pv_parks[0].surf_pv)*dt/e_ratio
            
            # Hydro
            p_m['dam'][i] = sum(cell.dam.p_out[m[i-1]:m[i]])*dt/e_ratio
            p_m['runriver'][i] = sum(cell.runriver.p_gen[m[i-1]:m[i]])*dt/e_ratio
            # p_m['phes_out'][i] = sum(cell.phes.p_out[m[i-1]:m[i]])*dt/e_ratio
            # p_m['phes_in'][i] = sum(cell.phes.p_in[m[i-1]:m[i]])*dt/e_ratio
            
            # Conventional
            p_m['nuc'][i] = sum(cell.nuclear.p_gen[m[i-1]:m[i]])*dt/e_ratio
            p_m['coal'][i] = sum(cell.coal.p_gen[m[i-1]:m[i]])*dt/e_ratio
            # p_m['ccgt'][i] = sum(cell.ccgt.p_gen[m[i-1]:m[i]])*dt/e_ratio
            
            # Batteries
            # p_m['bat_out'][i] = sum(cell.bat_parks[0].p_out[m[i-1]:m[i]])*dt/e_ratio
            # p_m['bat_in'][i] = sum(cell.bat_parks[0].p_in[m[i-1]:m[i]])*dt/e_ratio
            
            # PtG
            # p_m['ptg'][i] = sum(cell.ptg.p_in[m[i-1]:m[i]])*dt/e_ratio
            
            # Exchange 
            #TODO: to change, calculer import_export global en chaque temps t, tout ce qui est 
            # positif est import tout ce qui est negatif est export, puis sommer par mois
            imp = 0
            exp = 0
            for neighbour in cells:
                imp += grid.e_exch[neighbour.name][cell.name][m[i-1]:m[i]]*dt/e_ratio
                exp += grid.e_exch[cell.name][neighbour.name][m[i-1]:m[i]]*dt/e_ratio
                
            p_m['e_exch_out'][i] = sum(np.minimum(exp - imp,0))
            p_m['e_exch_in'][i] = sum(np.maximum(imp-exp,0))
            
            # Import
            p_m['e_imp'][i] = sum(grid.e_import[cell.name][m[i-1]:m[i]])*dt/e_ratio
            
            # Consumption
            p_m['cons'][i] = sum(cell.cons[m[i-1]:m[i]])*dt/e_ratio
                    
                    
        p_m_cells[cell.name] = p_m
        
        
        
        """ Plots """
        figure(fig_n, figsize = (7,4))
        
        
        N = 12
        ind = np.arange(N)    # the x locations for the groups
        width = 0.35       # the width of the bars: can also be len(x) sequence
        
        # Storage and export
        btm = 0
        # p0 = bar(ind, p_m['bat_in'], width, bottom = btm, color = '#bfbfbf', label = 'Battery in') #grey
        # btm += p_m['bat_in']
        # p1 = bar(ind, p_m['ptg'], width, bottom = btm, color = '#cc3300', label = 'PtG') #red
        # btm += p_m['ptg']
        # p2 = bar(ind, p_m['phes_in'], width, bottom = btm, color = '#2791f4', label = 'PHES in') #light blue
        # btm += p_m['phes_in']
        # p3 = bar(ind, p_m['e_exch_out'], width, bottom = btm, color = '#000000', label = 'Export') #black
        
        # Production
        btm = 0
        p4 = bar(ind, p_m['pv'], width, bottom = btm, color = '#ffe402', label = 'PV') # yellow
        btm += p_m['pv']
        p5 = bar(ind, p_m['wt_on'], width, bottom=btm, color = '#a3e2a1', label = 'WT onshore')
        btm += p_m['wt_on']
        p6 = bar(ind, p_m['wt_off'], width, bottom=btm, color = '#00B54E', label = 'WT offshore')
        btm += p_m['wt_off']
        p7 = bar(ind, p_m['dam'], width, bottom=btm, color = '#0003c7', label = 'Dam')
        btm += p_m['dam']
        p8 = bar(ind, p_m['runriver'], width, bottom=btm, color = '#42FFFC', label = 'Runriver')
        btm +=p_m['runriver']
        # p10 = bar(ind, p_m['nuc'], width, bottom=btm, color = '#77ffc8', label = 'Nuclear')
        # btm += p_m['nuc']
        # p11 = bar(ind, p_m['coal'], width, bottom=btm, color = '#7A3B00', label = 'Coal')
        # btm += p_m['coal']
        # p13 = bar(ind, p_m['bat_out'], width, bottom=btm, color = '#606060', label = 'Battery out')
        # btm += p_m['bat_out']
        # p12 = bar(ind, p_m['ccgt'], width, bottom=btm, color = '#ff7373', label = 'CCGT')
        # btm += p_m['ccgt']
        # p9 = bar(ind, p_m['phes_out'], width, bottom=btm , color = '#8ce2ff' , label = 'PHES out'')
        # btm += p_m['phes_out']
        p14 = bar(ind, p_m['e_exch_in'], width, bottom=btm , color = '#000000', label = 'Import')
        btm += p_m['e_exch_in']
        # p15 = bar(ind, p_m['e_imp'], width, bottom=btm, color = '#FFBB87', label = 'Import (out system)')
        
        #Consumption
        
        plot(ind, p_m['cons'], '--o', color='r', label ='Consumption')
        
        
        ylabel('[TWh]')
        title('Monthly energy flows in '+ name[cell.name])
        xticks(ind, ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
        #yticks(np.arange(0, 81, 10))
        
        legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=3)
        
        subplots_adjust(bottom=0.27, top = 0.91)
        
        savefig('img/'+ cell.name+'_prod_monthly.eps', format='eps', dpi=1000)
        
        fig_n += 1
        
    return fig_n