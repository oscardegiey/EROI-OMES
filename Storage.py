import numpy as np

class Battery:
    def __init__(self, specs, name, size, dt):
        self.name = name
        self.size = size
        self.dt = dt
        
        self.Esizemin = specs['Esizemin']        
        self.Esizemax = specs['Esizemax']
        self.efficiency_in = specs['efficiency_in']
        self.efficiency_out = specs['efficiency_out']
        self.alpha_in = specs['alpha_in']
        self.alpha_out = specs['alpha_out']
        self.leak = specs['leak'] / 24
        self.DoD = specs['DoD']
        self.n_cycle = specs['n_cycle']
        self.lifetime = specs['lifetime']*8760
        self.eps_inst = specs['eps_inst']
        self.eps_ut = - self.eps_inst * self.efficiency_in * self.dt / (self.DoD * self.n_cycle)
        
    def set_opti(self, capacity, energy, p_in, p_out):
        self.capacity = capacity            # [GWh] Optimal storage capacity
        self.energy = energy                # [GWh] Energy stored
        self.p_in = p_in
        self.p_out = p_out
        
class PHES:
    def __init__(self, specs, size):
        self.size = size
        self.p_inst_in = specs['p_inst_in']
        self.p_inst_out = specs['p_inst_out']
        self.capacity = specs['capacity']
        self.CO2perMWh = specs['CO2perMWh']
        self.local_in = np.zeros(size)
        self.efficiency_out = specs['efficiency_out']
        self.efficiency_in = specs['efficiency_in']
        self.ESOI = specs['esoi']

    def set_opti(self, energy, p_in, p_out): #Optimised prod is assigned
        self.p_out=p_out #output power
        self.p_in = p_in  #input power
        #self.p_pump = p_pump
        self.energy = energy #variation of the energy in the dam
        
class PtG:
    def __init__(self, specs, size, dt):
        self.size = size
        self.dt = dt
        self.Pmin = specs['p_min']
        self.Pmax = specs['p_max']
        self.efficiency_in = specs['efficiency_in']
        self.n_cycle = specs['n_cycle']
        self.eps_inst = specs['eps_inst']
        self.alpha = specs['alpha']
        self.eps_ut = - self.eps_inst * self.efficiency_in * self.dt / self.n_cycle
        self.lifetime = specs['lifetime']*8760
                
    def set_opti(self, p_inst, p_in):
        self.p_inst = p_inst           # Optimal installed power
        self.p_in = p_in               # Power input
        
class Gas:
    def __init__(self, specs, size):
        self.size = size
        self.Esizemin = specs['Esizemin']
        self.Esizemax = specs['Esizemax']
        self.eps = specs['eps']
        self.lifetime = specs['lifetime']*8760
        
    def set_opti(self, capacity, energy):
        self.capacity = capacity            # Optimal storage capacity
        self.energy = energy                # Energy stored
        
