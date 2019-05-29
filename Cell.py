try:
    import ruamel_yaml as yaml # Installed via Conda
except:
    try:
        import ruamel.yaml as yaml # Installed via pip
    except:
        print("ruamel.yaml not installed")
import numpy as np
from Producers import *
from Storage import *
import sys
from pulp import *
import matplotlib.pyplot as plt
from colors import colors as clr

class Cell:
    """Cell generated from a cell configuration file

    Attributes
    ----------
    size : int
        The size of the problem
    dt : int
        The timestep of the simulation, in hours


    """
    def __init__(self, size, dt, log_level):
        self.size      = size
        self.dt        = dt
        self.log_level = log_level

    def configure(self, config_file):
        stream = open(config_file, 'r')
        config = yaml.safe_load(stream)

        self.name      = config['name']
        self.re_share  = config['re_share']
        self.imp_share = config['re_share']
        self.nuc_share = config['nuc_share']
        if self.log_level >= 1:
            print("Loading cell: {0}".format(self.name))

        # Sets an order of importance for typical days to approximate the different inputs
        self.weight= config['weight']

        # Set consumption
        if self.log_level == 2: print('\tSet consumption,', end='')
        self.cons = np.loadtxt(open(config['consumption_file'],'r'), skiprows=1)
        self.cons = self.cons.reshape(-1,self.dt).mean(axis=1) # Reshape according to the time step
        self.cons = self.cons * 1e6

        # Producers
        self.wt_parks = []
        self.pv_parks = []

        ## PV
        if self.log_level == 2: print('Set PV parks,', end='')
        if 'regions' in config['pv']:
            for i in range(len(config['pv']['regions'])):
                pvpark = PVPark(config['pv']['regions'][i], config['pv'], self.size, self.dt)
                self.pv_parks.append(pvpark)
        else:
            pvpark = PVPark(config['pv'], config['pv'], self.size, self.dt)
            self.pv_parks.append(pvpark)

        ## WT
        if self.log_level == 2: print('Set WT parks,', end='')
        for onoff in config['wind']:
            if 'regions' in config['wind'][onoff]:
                for i in range(len(config['wind'][onoff]['regions'])):
                    wtpark = WTPark(config['wind'][onoff]['regions'][i], config['wind'][onoff], self.size, self.dt, onoff == 'onshore')
                    self.wt_parks.append(wtpark)
            else:
                wtpark = WTPark(config['wind'][onoff], config['wind'][onoff], self.size, self.dt, onoff == 'onshore')
                self.wt_parks.append(wtpark)

        ## CCGT
        if self.log_level == 2: print('Set CCGT,', end='')
        self.ccgt = CCGTPlant(config['ccgt'])

        ## Coal
        if self.log_level == 2: print('Set Coal,', end='')
        self.coal = CoalPlant(config['coal'])

        ## Nuclear
        if self.log_level == 2: print('Set Nuclear,', end='')
        self.nuclear = NuclearPlant(config['nuclear'], self.size)

        ## Hydro Dam
        if self.log_level == 2: print('Set Dam,', end='')
        self.dam = DamPlant(config['hydro_dam'], self.size, self.dt)

        ## Run-of-the-river
        if self.log_level == 2: print('Set Run-of-the-river,', end='')
        self.runriver = RunRiver(config['hydro_river'], self.size, self.dt)


        # Storage
        self.bat_parks = []

        ## Batteries
        if self.log_level == 2: print('Set Battery parks,', end='')
        batt = Battery(config['battery'], "Lithium", self.size, self.dt)
        self.bat_parks.append(batt)

        ## PHES
        if self.log_level == 2: print('Set PHES,', end='')
        self.phes = PHES(config['phes'], self.size)

        ## PtG
        if self.log_level == 2: print('Set Power-to-Gas,', end='')
        self.ptg = PtG(config['ptg'], self.size, self.dt)

        ## Gas
        if self.log_level == 2: print('Set Gas storage')
        self.gas = Gas(config['gas'], self.size)

    def set_opti(self, v, e_exch, gas_import, case):
        dtperday = int(24/case.dt)
        if case.optiTD:
            TDarray = [None] * case.size
            simpleTD = [None] * case.nrTD
            runner = 0
            hTD = []
            TDays = case.TDays
            for p in range(case.size):
                if TDays[p] not in simpleTD:
                    simpleTD[runner] = TDays[p]
                    runner += 1
            simpleTD.sort()

            for p in range(case.size):
                TDarray[p] = simpleTD.index(TDays[p])
            for t in range(365):
                hTD[t*dtperday:t*dtperday+dtperday] = range(dtperday)

        # Cell consumption to TD consumption
        cons = np.zeros(self.size)
        if case.optiTD:
            for t in range(self.size):
                h = int(hTD[t])
                TD = int(TDarray[t])
                cons[t] = self.cons[dtperday * (int(TDays[t])-1) + h]
            self.cons = cons

        for park in self.pv_parks:
            surf_pv = value(v[park.name + '_surf_pv'])
            p_gen = np.zeros(self.size)
            if case.optiTD:
                for t in range(self.size):
                    h = int(hTD[t])
                    TD = int(TDarray[t])
                    p_gen[t] = surf_pv * park.prod_unit[dtperday*(int(TDays[t])-1)+h]
            else:
                p_gen = surf_pv * park.prod_unit
            park.set_opti(surf_pv, p_gen)
            if self.log_level == 2:
                print('\tPV park {}, surface : {:.3f} km^2 [{:.3f} ; {:.3f}]'.format(park.name, surf_pv/1e6, park.minM2/1e6, park.maxM2/1e6))

        for park in self.wt_parks:
            n_wt = value(v[park.name + '_n_wt'])
            p_gen = np.zeros(self.size)
            if case.optiTD:
                for t in range(self.size):
                    h = int(hTD[t])
                    TD = int(TDarray[t])
                    p_gen[t] = n_wt * park.prod_unit[dtperday*(int(TDays[t])-1)+h]
            else:
                p_gen = n_wt * park.prod_unit
            park.set_opti(n_wt, p_gen)
            if self.log_level == 2:
                print('\tWT park {}, number of WT : {:.3f} [{:.3f} ; {:.3f}]'.format(park.name, n_wt, park.n_min, park.n_max))

        ccgt = np.zeros(self.size)
        coal = np.zeros(self.size)
        p_out_dam = np.zeros(self.size)
        e_dam = np.zeros(self.size)
        p_gen_runriver = np.zeros(self.size)

        for t in range(self.size):
            if case.optiTD:
                h = int(hTD[t])
                TD = int(TDarray[t])
                ccgt[t] = value(v['p_ccgt'][h+dtperday*TD])
                coal[t] = value(v['p_coal'][h+dtperday*TD])
                p_out_dam[t] = value(v['p_out_dam'][h+dtperday*TD])
                p_gen_runriver[t] = value(v['p_gen_runriver'][h+dtperday*TD])
            else:
                ccgt[t] = value(v['p_ccgt'][t])
                coal[t] = value(v['p_coal'][t])
                p_out_dam[t] = value(v['p_out_dam'][t])
                p_gen_runriver[t] = value(v['p_gen_runriver'][t])
            e_dam[t] = value(v['e_dam'][t])

        self.ccgt.set_opti(ccgt, gas_import)
        self.coal.set_opti(coal)
        self.dam.set_opti(p_out_dam,e_dam)
        self.runriver.set_opti(p_gen_runriver)

        nuclear = value(v['p_nuc'])
        self.nuclear.set_opti(nuclear)

        for park in self.bat_parks:
            capacity = value(v[park.name+'_capacity_bat'])
            energy = np.zeros(self.size)
            p_in = np.zeros(self.size)
            p_out = np.zeros(self.size)
            for t in range(self.size):
                energy[t] = value(v[park.name+'_e_bat'][t])
                if case.optiTD==1:
                    h = int(hTD[t])
                    TD = int(TDarray[t])
                    p_in[t] = value(v[park.name+'_p_in_bat'][h+dtperday*TD])
                    p_out[t] = value(v[park.name+'_p_out_bat'][h+dtperday*TD])
                else:
                    p_in[t] = value(v[park.name+'_p_in_bat'][t])
                    p_out[t] = value(v[park.name+'_p_out_bat'][t])

            park.set_opti(capacity,energy,p_in,p_out)
            if self.log_level == 2:
                print("\tBattery park {}, capacity : {:.3f} GWh".format(park.name, capacity/1e9))

        # PHES
        energy = np.zeros(self.size)
        p_in = np.zeros(self.size)
        p_out = np.zeros(self.size)
        for t in range(self.size):
            if case.optiTD==1:
                h = int(hTD[t])
                TD = int(TDarray[t])
                p_in[t] = value(v['p_in_phes'][h+dtperday*TD])
                p_out[t] = value(v['p_out_phes'][h+dtperday*TD])

            else:
                p_in[t] = value(v['p_in_phes'][t])
                p_out[t] = value(v['p_out_phes'][t])
            energy[t] = value(v['e_phes'][t])

        self.phes.set_opti(energy,p_in,p_out)

        # PtG
        p_inst = value(v['p_inst_ptg'])
        p_in = np.zeros(self.size)
        for t in range(self.size):
            if case.optiTD==1:
                h = int(hTD[t])
                TD = int(TDarray[t])
                p_in[t] = value(v['p_in_ptg'][h+dtperday*TD])

            else:
                p_in[t] = value(v['p_in_ptg'][t])

        self.ptg.set_opti(p_inst,p_in)
        if self.log_level == 2:
            print("\tPower-to-Gas, installed power : {:.3f} GW".format(p_inst/1e9))

        capacity = value(v['capacity_gas'])
        energy = np.zeros(self.size)
        for t in range(self.size):
            energy[t] = value(v['e_gas'][t])
        self.gas.set_opti(capacity,energy)
        if self.log_level == 2:
            print("\tGas storage, capacity : {:.3f} TWh".format(capacity/1e12))

        self.e_exch = e_exch[self.name]

    def display(self, disp, timeframe=[0,8760], resolution="monthly", ax=None):
    # timeframe is expressed in [start,end] with start and end between [0,8760]
        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
        x = np.arange(0, self.size)
        if disp == "simple_balance":
            PV = sum(park.p_gen for park in self.pv_parks)
            WT = sum(park.p_gen for park in self.wt_parks)
            Dam = self.dam.p_out
            Runriver = self.runriver.p_gen
            RE = (PV + WT + Dam + Runriver) / 1e9

            Coal = self.coal.p_gen
            Nuclear = self.nuclear.p_gen
            CCGT = self.ccgt.p_NRE
            NRE = (Coal + Nuclear + CCGT) / 1e9

            Bat = sum(park.p_out for park in self.bat_parks)
            PHES = self.phes.p_out
            PtG = self.ccgt.p_RE
            PROD_STORAGE = (Bat + PHES + PtG) / 1e9

            Bat = sum(park.p_in for park in self.bat_parks)
            PHES = self.phes.p_in
            PtG = self.ptg.p_in
            CONS_STORAGE = (Bat + PHES + PtG) / 1e9
            EXCH = sum(self.e_exch[neigh] for neigh in self.e_exch) / 1e9

            CURTAILMENT = RE + NRE + PROD_STORAGE + CONS_STORAGE - self.cons / 1e9 - EXCH

            if np.sum(NRE) == 0:
                ax.stackplot(x, RE, PROD_STORAGE, -EXCH, labels=['RE','STOR', 'EXCH'])
            else:
                ax.stackplot(x, RE, NRE, PROD_STORAGE, -EXCH, labels=['RE','NRE','STOR', 'EXCH'])
            ax.plot(x,self.cons / 1e9, label='Consumption')
            ax.plot(x, CURTAILMENT, label='Curtailment')

            ax.set_xlim(timeframe[0]/self.dt,timeframe[1]/self.dt)
            ax.legend()
            ax.grid(True)
            ax.set_ylabel('Power in GW')
            ax.set_xlabel('Time')

        elif disp == "complex_balance":
            PV = sum(park.p_gen for park in self.pv_parks)
            WTon = np.zeros(self.size); WToff = np.zeros(self.size)
            for park in self.wt_parks:
                if park.onshore:
                    WTon += park.p_gen
                else:
                    WToff += park.p_gen
            Dam = self.dam.p_out
            Runriver = self.runriver.p_gen
            RE = (PV + WTon + WToff + Dam + Runriver)

            Bat = sum(park.p_out for park in self.bat_parks)
            PHES = self.phes.p_out
            PtG = self.ccgt.p_RE
            PROD_STORAGE = (Bat + PHES + PtG)

            EXCH = np.zeros(self.size)
            for neighbour in self.e_exch:
                EXCH += -self.e_exch[neighbour]
            EXCH = np.maximum(0,EXCH)

            ax.stackplot(x, PV, WTon, WToff, Dam, Runriver, Bat, PHES, PtG, EXCH, labels=['PV','WTon','WToff','Dam','Runriver','Bat','PHES','PtG','Exchange'], colors=(clr('PV'),clr('WTon'),clr('WToff'),clr('Dam'),clr('Runriver'),clr('Battery'),clr('PHES'),clr('PtG'),clr('Exchange')))
            ax.plot(x,RE, '--g',label='RE')
            ax.plot(x,self.cons, '-r',linewidth=2,label='Consumption')
            ax.set_xlim(timeframe[0]/self.dt,timeframe[1]/self.dt)
            ax.set_ylim(0,max(self.cons))
            ax.legend()
            ax.grid(True)
            ax.set_ylabel('Power in W')
            ax.set_xlabel('Time')
            ax.set_title('Energy balance of ' + self.name)

        elif disp == "flows":
            if resolution == "dt":
                dtin = np.arange(timeframe[0] / self.dt, timeframe[1] / self.dt + 1).astype(int)
                x = np.arange(timeframe[0] / self.dt, timeframe[1] / self.dt)
            elif resolution == "daily":
                dtin = np.arange(timeframe[0] / self.dt, (timeframe[1]+24) / self.dt, 24 / self.dt).astype(int)
                x = np.arange(timeframe[0] / self.dt, timeframe[1] / self.dt, 24 / self.dt) / 24
            elif resolution == "weekly":
                dtin = np.arange(timeframe[0] / self.dt, (timeframe[1]+168) / self.dt, 168 / self.dt).astype(int)
                x = np.arange(timeframe[0] / self.dt, timeframe[1] / self.dt, 168 / self.dt) / 168
            elif resolution == "monthly":
                din = np.array([0,31,28,31,30,31,30,31,31,30,31,30,31])
                dtin = np.cumsum(din * 24 / self.dt).astype(int)
                x = np.arange(12)

            size = len(dtin)-1

            PV = np.zeros(size); WTon = np.zeros(size); WToff = np.zeros(size)
            Dam = np.zeros(size); Runriver = np.zeros(size); Bat = np.zeros(size)
            Batout = np.zeros(size); Batin = np.zeros(size); CCGT = np.zeros(size)
            PtG = np.zeros(size); PHES = np.zeros(size); PHESout = np.zeros(size)
            PHESin = np.zeros(size); Exchange = np.zeros(size); Exchangeout = np.zeros(size)
            Exchangein = np.zeros(size); Cons = np.zeros(size)

            for i in range(len(dtin)-1):
                start = dtin[i]
                end = dtin[i+1]
                PV[i] = sum(sum(park.p_gen[start:end]) * self.dt for park in self.pv_parks)
                for park in self.wt_parks:
                    if park.onshore:
                        WTon[i] = sum(park.p_gen[start:end]) * self.dt
                    else:
                        WToff[i] = sum(park.p_gen[start:end]) * self.dt
                Dam[i] = sum(self.dam.p_out[start:end]) * self.dt
                Runriver[i] = sum(self.runriver.p_gen[start:end]) * self.dt
                Bat[i] = sum(sum(park.p_out[start:end]) + sum(park.p_in[start:end]) * self.dt for park in self.bat_parks)
                CCGT[i] = sum(self.ccgt.p_RE[start:end]) * self.dt
                PtG[i] = sum(self.ptg.p_in[start:end]) * self.dt
                PHES[i] = (sum(self.phes.p_out[start:end]) + sum(self.phes.p_in[start:end])) * self.dt
                for neighbour in self.e_exch:
                    Exchange[i] += -sum(self.e_exch[neighbour][start:end]) * self.dt
                Cons[i] = sum(self.cons[start:end]) * self.dt

            Batout = np.maximum(0,Bat)
            Batin = np.minimum(0,Bat)
            PHESout = np.maximum(0,PHES)
            PHESin = np.minimum(0,PHES)
            Exchangeout = np.maximum(0,Exchange)
            Exchangein = np.minimum(0,Exchange)

            ax.bar(x,PV, color=clr('PV'), label='PV'); btm = PV
            ax.bar(x,WTon, color=clr('WTon'), label='WTon', bottom=btm); btm += WTon
            ax.bar(x,WToff, color=clr('WToff'), label='WToff', bottom=btm); btm += WToff
            ax.bar(x,Dam, color=clr('Dam'), label='Dam', bottom=btm); btm += Dam
            ax.bar(x,Runriver, color=clr('Runriver'), label='Runriver', bottom=btm); btm += Runriver
            ax.bar(x,Batout, color=clr('Battery'), label='Battery', bottom=btm); btm += Batout
            ax.bar(x,CCGT, color=clr('CCGT'), label='CCGT', bottom=btm); btm += CCGT
            ax.bar(x,PHESout, color=clr('PHES'), label='PHES', bottom=btm); btm += PHESout
            ax.bar(x,Exchangeout, color=clr('Exchange'), label='Exchange', bottom=btm)

            ax.bar(x,Batin, color=clr('Battery')); btm = Batin
            ax.bar(x,PtG, color=clr('PtG'), label='PtG', bottom=btm); btm += PtG
            ax.bar(x,PHESin, color=clr('PHES'), bottom=btm); btm += PHESin
            ax.bar(x,Exchangein, color=clr('Exchange'), bottom=btm)

            ax.plot(x,Cons, '--o', color='red', label='Consumption')

            ax.set_ylabel('Energy in Wh')
            ax.set_xlabel('Weeks')
            #ax.set_xlim(timeframe[0]/(8760/size)-0.5,timeframe[1]/(8760/size)-0.5)
            ax.axhline(0, linestyle='--', color='black')
            ax.legend()
            ax.set_title('Weekly energy flows of ' + self.name)

        elif disp == 'storages':
            x = np.arange(self.size)
            PV = sum(park.p_gen for park in self.pv_parks)
            WT = sum(park.p_gen for park in self.wt_parks)
            Dam = self.dam.p_out
            Runriver = self.runriver.p_gen
            RE = (PV + WT + Dam + Runriver) / 1e9

            Coal = self.coal.p_gen
            Nuclear = self.nuclear.p_gen
            CCGT = self.ccgt.p_NRE
            NRE = (Coal + Nuclear + CCGT) / 1e9

            Bat = sum(park.p_out for park in self.bat_parks)
            PHES = self.phes.p_out
            PtG = self.ccgt.p_RE
            PROD_STORAGE = (Bat + PHES + PtG) / 1e9

            Bat = sum(park.p_in for park in self.bat_parks)
            PHES = self.phes.p_in
            PtG = self.ptg.p_in
            CONS_STORAGE = (Bat + PHES + PtG) / 1e9

            EXCH = sum(self.e_exch[neigh] for neigh in self.e_exch) / 1e9

            curt = RE + NRE + PROD_STORAGE + CONS_STORAGE - self.cons / 1e9 - EXCH
            CURTAILMENT = self.cons / 1e9 + EXCH - CONS_STORAGE - RE - NRE - PROD_STORAGE
            ax.plot(x, curt)
            #ax.plot(x, self.bat_parks[0].energy/self.bat_parks[0].capacity)
            #ax.plot(x, self.gas.energy / self.gas.capacity)
            ax.set_xlim(timeframe[0]/self.dt,timeframe[1]/self.dt)
            ax.grid(True)


        return ax
