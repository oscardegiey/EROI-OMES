try:
    import ruamel_yaml as yaml # Installed via Conda
except:
    try:
        import ruamel.yaml as yaml # Installed via pip
    except:
        print("ruamel.yaml not installed")

from Cell import Cell
from Grid import Grid
from pulp import *
from KMedoids import KMedoids
import numpy as np
import os, glob, multiprocessing, time, json
import matplotlib.pyplot as plt
from colorama import Style

def serialize(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj.__dict__
    return obj

class Case:
    """Case generated from a simulation configuration file.

    Attributes
    ----------
    log_level : int
        A value between 0 and 2: 0 for little information, 2 for lots of
        information.
    name : str
        The name of the case
    dt : int
        The timestep of the simulation, in hours
    size : int
        The size of the problem
    use_CPLEX : bool
        1 if CPLEX solver is used, 0 if CBC solver is used
    optiTD : bool
        1 if Typical Days are used, 0 if not
    replicate : int
        The number of times the K-Medoids algorithm is used to find the best
        Typical Days
    loadTD : bool
        1 if Typical Days are loaded instead of computed, 0 if not
    nrTD : int
        The number of Typical Days used
    e_exch : array
        2D array representing the maximum electricity exchange between two
        cells in W
    gas_exch : array
        2D array representing the maximum gas exchange between two cells in W
    e_imp : array
        1D array representing the maximum electricity import of each cell in W
    gas_imp : array
        1D array representing the maximum gas import of each cell in W
    grid : object
        Grid object corresponding to the case
    TDays : array
        1D array representing the corresponding typical day for every timestep
    """

    def __init__(self, configuration_file, log_level):
        """Create Case object based on the given configuration file.

        Parameters
        ----------
        configuration_file: str
            Path to the case configuration file (a .yml file)
        log_level: int
            A value between 0 and 2: 0 for little information, 2 for lots of
            information

        Returns
        -------
        object
            Case object

        """

        self.log_level = log_level
        stream         = open(configuration_file, 'r')
        config         = yaml.safe_load(stream)



        # === Load parameters ===
        self.name      = config['name']
        self.dt        = config['dt']
        self.size      = int(8760/self.dt)
        self.use_CPLEX = config['use_CPLEX']
        self.optiTD    = config['optiTD']

        if self.log_level >= 1:
            print(Style.BRIGHT + "=== Case configuration ===" + Style.NORMAL)
        if self.log_level == 2:
            print("Case name: \t{0}\nTime step: \t{1}\nProblem size: \t{2}\nCPLEX solver: \t{3}\nTypical Days: \t{4}".format(self.name, self.dt, self.size, self.use_CPLEX, self.optiTD))

        if self.optiTD: # If Typical Days are used
            self.replicate = config['replicate']
            self.loadTD    = config['loadTD']
            self.nrTD      = config['nrTD']
            if self.log_level == 2:
                print("Load TD: \t{0}\nNumber of TD: \t{1}".format(self.loadTD, self.nrTD))
                if not(self.loadTD):
                    print("Replicate: \t{0}".format(self.replicate))

        # Policies
        re_share      = config['re_share']
        imp_share     = config['imp_share']
        self.policies = [re_share, imp_share]
        if self.log_level == 2:
            print("Max RE share: \t{0}\nMax imp share: \t{1}".format(re_share,imp_share))

        # Cells
        self.cells = []
        nbcells = len(config['cells'])
        if self.log_level >= 1:
            print('\n' + Style.BRIGHT + "=== Cells configuration ===" + Style.NORMAL)
        for cell in config['cells']:
             self.cells.append(Cell(self.size, self.dt, self.log_level))
             self.cells[-1].configure(cell)

        # e_exch, gas_exch, e_imp and gas_imp
        e_exch = {}; gas_exch = {}; e_imp    = {}; gas_imp  = {}
        for i, cell in enumerate(self.cells):
            e_exch[cell.name]   = {}
            gas_exch[cell.name] = {}
            e_imp[cell.name]    = config['e_imp'][i]
            gas_imp[cell.name]  = config['gas_imp'][i]
            for j, neighbour in enumerate(self.cells):
                e_exch[cell.name][neighbour.name]  = config['e_exch'][i][j]
                gas_exch[cell.name][neighbour.name] = config['gas_exch'][i][j]

        self.e_exch   = e_exch
        self.gas_exch = gas_exch
        self.e_imp    = e_imp
        self.gas_imp  = gas_imp

        self.grid = Grid(self.e_imp, self.e_exch, self.gas_imp, self.gas_exch)

        # Computed values
        self.E            = {}; self.E_tot        = {}
        self.EROI_pv      = {}; self.EROI_wt      = {}
        self.share        = {}; self.share_tot    = {}
        self.LF           = {}; self.E_invest     = {}
        self.E_invest_tot = {};

        # Typical Days
        if self.optiTD:
            assert self.dt < 12, "Timestep to big to compute typical days."

            if self.loadTD: # Load Typical Days from least error file
                filename = os.path.join('typical_days', self.name + '_' +  \
                        'TDays'+ '_' + str(self.dt) + '_' + str(self.nrTD) \
                        + 'TD_*')
                possible_TD = glob.glob(filename)
                min_error = 1e9
                for f in possible_TD:
                    error = float(f.split('_')[-1][:-4])
                    if error < min_error:
                        TDays_f = f; min_error = error
                if self.log_level >= 1:
                    print("\n" + Style.BRIGHT + "=== Typical Days ===" + Style.NORMAL)
                    print('Min error : \t{:.6f}'.format(min_error))
                self.TDays = np.loadtxt(TDays_f)

            else: # Start multithreaded computation to find the Typical Days
                threads = 4
                if self.log_level >= 1:
                    print("\n" + Style.BRIGHT + "=== Typical Days ===" + Style.NORMAL)
                    print("Compute typical days with "+str(threads)+" threads")
                manager = multiprocessing.Manager()
                return_dict = manager.dict()
                jobs = [None] * threads

                self.replicate = int(self.replicate/threads)
                for i in range(threads): # Start threads
                    jobs[i] = multiprocessing.Process(target=self.computeTD, \
                                                      args=(i, return_dict))
                    jobs[i].start()
                for i in range(threads): # Wait for threads to finish
                    jobs[i].join()

                best_error = 1e6
                for (TDays, error) in return_dict.values():
                    if error < best_error:
                        self.TDays = TDays.astype('int')
                        best_error = error
                if self.log_level >= 1:
                    print('Min error : \t{:.6f}'.format(best_error))

                dtperday = int(24 / self.dt)
                dtpd = int(24 / self.dt)
                powerTD = np.zeros(self.size)

                # Plot the load duration curve with and without typical days
                # for the user to assess the precision
                if self.log_level == 2:
                    for cell in self.cells:
                        for i in range(0, self.size, dtperday):
                            TD = self.TDays[i]
                            powerTD[i:i+dtpd] = cell.cons[TD*dtpd:TD*dtpd+dtpd]

                        sort = np.sort(cell.cons)[::-1]
                        exceedence = np.arange(1.,len(sort)+1) / len(sort)
                        l1, = plt.plot(exceedence*100, sort)
                        sort = np.sort(powerTD)[::-1]
                        exceedence = np.arange(1.,len(sort)+1) / len(sort)
                        l2, = plt.plot(exceedence*100, sort)
                        plt.legend((l1,l2), ('Without TD', 'With TD'))
                        plt.title("LDC with and without Typical Days : " + cell.name)
                        plt.show()


                filename = os.path.join('typical_days', self.name + '_' +    \
                        'TDays'+ '_' + str(self.dt) + '_' + str(self.nrTD) + \
                        'TD_' + str(best_error) + '.txt')
                np.savetxt(filename, self.TDays, fmt='%d')

    def computeTD(self, procnum, results):
        """Compute Typical Days with K-Medoid Algorithm

        Parameters
        ----------
        procnum : int
            The number of the thread running this function
        results : list
            A list that contains the best Typical Days and their errors for
            each thread
        """
        dim = 2*len(self.cells) + sum(len(cell.pv_parks) + len(cell.wt_parks) \
                for cell in self.cells)
        distrib = np.zeros((self.size,dim))
        count   = 0

        # Construction and normalization of the matrix needed for the clustering
        for cell in self.cells:
            w = cell.weight
            distrib[:,count] = w[0] * cell.cons / np.sum(cell.cons); count += 1
            distrib[:,count] = w[1] * cell.runriver.local_in / np.sum(cell.runriver.local_in); count += 1
            for park in cell.pv_parks:
                distrib[:,count] = w[2] * park.irr / np.sum(park.irr); count += 1
            for park in cell.wt_parks:
                distrib[:,count] = w[3] * park.wind / sum(park.wind); count += 1

        KMedoidsMatrix = distrib.reshape((365,int(dim*24/self.dt))).tolist()
        best_medoids   = []
        best_clusters  = {}
        best_error     = 1e9
        tic            = time.time()

        for i in range(self.replicate):
            k_medoids = KMedoids(n_cluster=self.nrTD,max_iter=100,tol=0.0001)
            k_medoids.fit(KMedoidsMatrix)
            error = np.sum(list(k_medoids.cluster_distances.values()))/self.nrTD
            if error < best_error:
                best_medoids  = k_medoids.medoids
                best_clusters = k_medoids.clusters
                best_error    = error
            if (self.log_level == 2):
                if (i % max(1,int(self.replicate/100)) == 0 and procnum == 0):
                    time_remaining = (time.time() - tic) / (i+1) * (self.replicate - (i+1))
                    print('TD computation : ' + str(int(i/self.replicate*100)) + '%' + '\tTime remaining : {:10.3f} s'.format(time_remaining), end='\r')
                    #print('TD computation : ' + str(int(i/self.replicate*100)) + '%' + '\tTime remaining : {:10.3f} s'.format(time_remaining))
        if (self.log_level == 2 and procnum == 0):
            print('TD computation : 100%\tTotal time : {:10.3f} s             '.format(time.time() - tic))
        TDays        = np.zeros(self.size)
        dtperday     = int(24 / self.dt)
        best_medoids = [*best_medoids]
        best_medoids.sort()
        for med in best_medoids:
            for day in best_clusters[med]:
                TDays[dtperday*day:dtperday*day+dtperday] = med

        results[procnum] = (TDays, best_error)

    def set_opti(self,e_exch, e_import, gas_exch, gas_import, lpvariables):
        if self.log_level >= 1:
            print("\n" + Style.BRIGHT + "=== Set optimal solution ===" + Style.NORMAL)
        dtperday = int(24/self.dt)
        if self.optiTD:
            TDarray = [None] * self.size
            simpleTD = [None] * self.nrTD
            runner = 0
            hTD = []
            TDays = self.TDays
            for p in range(self.size):
                if TDays[p] not in simpleTD:
                    simpleTD[runner] = TDays[p]
                    runner += 1
            simpleTD.sort()

            for p in range(self.size):
                TDarray[p] = simpleTD.index(TDays[p])
            for t in range(365):
                hTD[t*dtperday:t*dtperday+dtperday] = range(dtperday)

        e_exch_opt = {}
        e_import_opt = {}
        gas_exch_opt = {}
        gas_import_opt = {}

        for cell in self.cells:
            e_exch_opt[cell.name] = {}
            gas_exch_opt[cell.name] = {}

            # -- Import --
            e_import_opt[cell.name] = np.zeros(self.size)
            gas_import_opt[cell.name] = np.zeros(self.size)
            for t in range(self.size):
                e_import_opt[cell.name][t] = value(e_import[cell.name][t])
                gas_import_opt[cell.name][t] = value(gas_import[cell.name][t])

            # -- Exchange --
            for neighbour in self.cells:
                e_exch_opt[cell.name][neighbour.name] = np.zeros(self.size)
                gas_exch_opt[cell.name][neighbour.name] = np.zeros(self.size)
                for t in range(self.size):
                    if self.optiTD:
                        h = int(hTD[t])
                        TD = int(TDarray[t])
                        e_exch_opt[cell.name][neighbour.name][t] = value(e_exch[cell.name][neighbour.name][h+dtperday*TD])
                        gas_exch_opt[cell.name][neighbour.name][t] = value(gas_exch[cell.name][neighbour.name][h+dtperday*TD])
                    else:
                        e_exch_opt[cell.name][neighbour.name][t] = value(e_exch[cell.name][neighbour.name][t])
                        gas_exch_opt[cell.name][neighbour.name][t] = value(gas_exch[cell.name][neighbour.name][t])

            v = lpvariables[cell.name]
            if self.log_level >= 1:
                print("Cell : {}".format(cell.name))
            cell.set_opti(v,e_exch_opt,gas_import_opt[cell.name], self)
            self.grid.set_opti(e_import_opt, e_exch_opt, gas_import_opt, gas_exch_opt)

    def compute_values(self):
        if self.log_level >= 1:
            print("\n" + Style.BRIGHT + "=== Compute values ===" + Style.NORMAL)
        self.E_tot['Prod'] = 0
        self.E_tot['RE'] = 0
        self.E_tot['Import'] = 0
        self.E_tot['Cons'] = 0

        for cell in self.cells:
            if self.log_level == 2:
                print("Cell : {}".format(cell.name))

            # -- Energy sums --
            self.E[cell.name] = {}
            self.E[cell.name]['Import'] = sum(self.grid.e_import[cell.name] * self.dt)
            self.E[cell.name]['Nuclear'] = sum(cell.nuclear.p_gen * self.dt)
            self.E[cell.name]['Coal'] = sum(cell.coal.p_gen * self.dt)
            self.E[cell.name]['CCGT'] = sum(cell.ccgt.p_gen * self.dt)
            self.E[cell.name]['Cons'] = sum(cell.cons * self.dt)
            self.E[cell.name]['PV'] = sum(sum(park.p_gen * self.dt) for park in cell.pv_parks)
            self.E[cell.name]['WT'] = sum(sum(park.p_gen * self.dt) for park in cell.wt_parks)
            self.E[cell.name]['Dam'] = sum(cell.dam.p_out * self.dt)
            self.E[cell.name]['Runriver'] = sum(cell.runriver.p_gen * self.dt)
            self.E[cell.name]['Gas_import'] = sum(self.grid.gas_import[cell.name]) * cell.ccgt.efficiency * self.dt

            self.E[cell.name]['RE'] = sum(self.E[cell.name][key] for key in ['Dam','Runriver','PV','WT'])
            self.E[cell.name]['Prod'] = sum(self.E[cell.name][key] for key in ['RE','Gas_import','Coal','Nuclear'])


            self.E[cell.name]['Batt_in'] = sum(sum(park.p_in for park in cell.bat_parks) * self.dt)
            self.E[cell.name]['Batt_out'] = sum(sum(park.p_out for park in cell.bat_parks) * self.dt)
            self.E[cell.name]['PHES_in'] = sum(cell.phes.p_in * self.dt)
            self.E[cell.name]['PHES_out'] = sum(cell.phes.p_out * self.dt)
            self.E[cell.name]['PtG_in'] = sum(cell.ptg.p_in * self.dt)
            self.E[cell.name]['Exch_in'] = sum(sum(np.maximum(0,self.grid.e_exch[neigh.name][cell.name]) * self.dt) for neigh in self.cells)
            self.E[cell.name]['Exch_out'] = sum(sum(np.maximum(0,self.grid.e_exch[cell.name][neigh.name]) *self.dt) for neigh in self.cells)

            self.E[cell.name]['Curtailment'] = sum(self.E[cell.name][key] for key in ['RE','CCGT','Coal','Nuclear','Batt_in','Batt_out','PHES_in','PHES_out','PtG_in','Exch_in']) \
                                             - sum(self.E[cell.name][key] for key in ['Exch_out','Cons'])


            # -- EROI PV --
            self.EROI_pv[cell.name] = {}

            # -- EROI PV --
            for park in cell.pv_parks:
                self.EROI_pv[cell.name][park.name] = sum(park.prod_unit * self.dt) / park.costpermsq * park.lifetime / (self.size * self.dt)

            # -- EROI WT --
            self.EROI_wt[cell.name] = {}
            for park in cell.wt_parks:
                self.EROI_wt[cell.name][park.name] = sum(park.prod_unit * self.dt) / park.costperWTinst * park.lifetime / (self.size * self.dt)

            # -- Shares --
            self.share[cell.name] = {}
            for key in ['RE','PV','WT','Runriver','Dam','Coal','Nuclear','CCGT','Import','Curtailment']:
                self.share[cell.name][key] = self.E[cell.name][key] / self.E[cell.name]['Prod']
            self.share[cell.name]['RE/cons'] = self.E[cell.name]['RE'] / self.E[cell.name]['Cons']
            if self.log_level == 2:
                print("\tRE share : \t{:.3f}\n\tCurtailment : \t{:.3f}".format(self.share[cell.name]['RE'], self.share[cell.name]['Curtailment']))

            # -- Load factors --
            self.LF[cell.name] = {}
            self.LF[cell.name]['Runriver'] = self.E[cell.name]['Runriver'] / (cell.runriver.p_inst * self.size * self.dt)
            self.LF[cell.name]['Dam'] = self.E[cell.name]['Dam'] / ((cell.dam.p_inst * self.size * self.dt) or 1)
            self.LF[cell.name]['PtG'] = -self.E[cell.name]['PtG_in'] / (cell.ptg.p_inst * self.size * self.dt)
            self.LF[cell.name]['Batt'] = -self.E[cell.name]['Batt_in'] / (sum(cell.bat_parks[i].capacity * cell.bat_parks[i].alpha_in * self.size * self.dt for i in range(len(cell.bat_parks))) or 1)
            self.LF[cell.name]['CCGT'] = self.E[cell.name]['CCGT'] / (cell.ccgt.Pmax * self.size * self.dt)

            # -- System Total --
            self.E_tot['Prod'] += self.E[cell.name]['Prod']
            self.E_tot['RE'] += self.E[cell.name]['RE']
            self.E_tot['Import'] += self.E[cell.name]['Import']
            self.E_tot['Cons'] += self.E[cell.name]['Cons']

        self.share_tot['RE'] = self.E_tot['RE'] / (self.E_tot['Prod'] + self.E_tot['Import'])

        self.E_invest['RE'] = {}
        self.E_invest['Storage'] = {}
        self.E_invest['RE0'] = {}
        self.E_invest_tot['RE'] = 0
        self.E_invest_tot['Storage'] = 0
        self.E_invest_tot['RE0'] = 0
        storage_weight = 0.5

        for cell in self.cells:
            dt = self.dt
            size = self.size

            # -- Energy invested in storage --
            self.E_invest['Storage'][cell.name] = 0
            self.E_invest['Storage'][cell.name] += -self.E[cell.name]['PHES_in'] / cell.phes.ESOI
            self.E_invest['Storage'][cell.name] += sum(cell.ptg.p_in)*cell.ptg.eps_ut*(1-storage_weight)
            self.E_invest['Storage'][cell.name] += cell.ptg.p_inst * cell.ptg.eps_inst * size * dt / cell.ptg.lifetime * storage_weight
            self.E_invest['Storage'][cell.name] += cell.gas.capacity * cell.gas.eps * size * dt / cell.gas.lifetime
            for park in cell.bat_parks:
                self.E_invest['Storage'][cell.name] += sum(park.p_in) * park.eps_ut * (1 - storage_weight)
                self.E_invest['Storage'][cell.name] += park.capacity * park.eps_inst * storage_weight * size * dt / park.lifetime

            # -- Energy invested in RE assets --
            self.E_invest['RE'][cell.name] = 0
            self.E_invest['RE0'][cell.name] = 0
            for park in cell.pv_parks:
                self.E_invest['RE'][cell.name] += park.surf_pv * park.costpermsq * size * dt / park.lifetime
                self.E_invest['RE0'][cell.name] += park.minM2 * park.costpermsq * size * dt / park.lifetime
            for park in cell.wt_parks:
                self.E_invest['RE'][cell.name] += park.n_wt * park.costperWTinst * size * dt / park.lifetime
                self.E_invest['RE0'][cell.name] += park.n_min * park.costperWTinst * size * dt / park.lifetime

            self.E_invest_tot['RE'] += self.E_invest['RE'][cell.name]
            self.E_invest_tot['Storage'] += self.E_invest['Storage'][cell.name]
            self.E_invest_tot['RE0'] += self.E_invest['RE0'][cell.name]

        self.EROI_tot = self.E_tot['Cons'] / (self.E_invest_tot['RE'] + self.E_invest_tot['Storage'])
        if self.log_level >= 1:
            print("System :")

            print("\tEROI : \t{:.3f}".format(self.EROI_tot))

    def save_json(self):
        a = json.dumps(self, default=serialize)

        f = open('CASERESULTS/' + str(self.optiTD) + 'TD_' + str(self.dt) + 'dt_REshare_' + str(self.policies[0]) + '.json','w')
        f.write(a)
        f.close()

        b = json.dumps(self.share, default=serialize)
        if self.optiTD:
            g = open('CASERESULTS/' + str(self.nrTD) + 'TD_' + str(self.dt) + 'dt_' +  'REshare_' + str(self.policies[0]) + '_Share' + '.json','w')
        else:
            g = open('CASERESULTS/' + str(self.optiTD) + 'TD_' + str(self.dt) + 'dt_' +  'REshare_' + str(self.policies[0]) + '_Share' + '.json','w')
        g.write(b)
        g.close()
