
import numpy as np


""" BATTERY """

class Battery:
    def __init__(self,type, size, dt):
        self.type = type # Type of battery used
        self.size = size
        self.dt = dt
        
    def set_spec(self, Esizemin,Esizemax, efficiency_in, efficiency_out, alpha_in, alpha_out, leak, DoD, n_cycle, lifetime, eps_inst):
        self.Esizemax=Esizemax              #[Wh] maximum storage size
        self.Esizemin=Esizemin              #[Wh] minimum storage size
        self.efficiency_in=efficiency_in    # efficiency of storage in
        self.efficiency_out=efficiency_out  # efficiency of storage out
        self.alpha_in = alpha_in   # charge rate [1/dt]
        self.alpha_out = alpha_out  # discharge rate [1/dt]
        self.leak=leak
        self.DoD=DoD                        # Deepth of discharge
        self.n_cycle=n_cycle                # Number of cycle which can be done
        self.lifetime=lifetime*8760         # Number of hours the storage can be used
        self.eps_inst=eps_inst              # Energy cost per energy unit : 100 <=> it cost 100 J to create 1 J of storage
        self.eps_ut = - eps_inst*efficiency_in*self.dt/(DoD*n_cycle) # Energy cost at utilisation
        
    def set_opti(self, capacity, energy, p_in, p_out):
        self.capacity = capacity            # [GWh] Optimal storage capacity
        self.energy = energy                # [GWh] Energy stored
        self.p_in = p_in
        self.p_out = p_out
        
""" PHES """

class PHES:
    def __init__(self, size):
        self.size = size
        
    def set_spec(self, p_inst_in, p_inst_out, capacity, CO2,local_in, eff_out, eff_in, ESOI):    
        self.p_inst_in = p_inst_in          # [W] Maximum input power
        self.p_inst_out = p_inst_out        # [W] Maximum output power
        self.capacity = capacity            # [Wh] Maximum energy stored
        self.CO2perMWh=CO2
        self.local_in = local_in            # [W] local inflow
        self.efficiency_out = eff_out       # efficiency of storage out
        self.efficiency_in = eff_in         # efficiency of storage in
        self.ESOI = ESOI                    # Storage ESOI
        
    def set_opti(self, energy, p_in, p_out): #Optimised prod is assigned
        self.p_out=p_out #output power
        self.p_in = p_in  #input power
        #self.p_pump = p_pump
        self.energy = energy #variation of the energy in the dam
        

""" PtG """


class PtG:
    def __init__(self,size, dt):
        self.size = size
        self.dt = dt
        
    def set_spec(self, Pmin, Pmax, efficiency_in, n_cycle, eps_inst, alpha, lifetime):
        self.Pmin = Pmin                    # PtG installed power
        self.Pmax= Pmax                     # PtG maximum powere
        self.efficiency_in=efficiency_in    # efficiency of storage in               
        self.n_cycle=n_cycle                # Number of cycle which can be done
        self.eps_inst=eps_inst              # [MWh/MW] Energy cost per energy unit : 100 <=> it cost 100 J to
        self.alpha = alpha                  # Minimal utilisation [ratio]
        self.eps_ut = - eps_inst*efficiency_in*self.dt/(n_cycle)
        self.lifetime = lifetime*8760
        
    def set_opti(self, p_inst, p_in):
        self.p_inst = p_inst           # Optimal installed power
        self.p_in = p_in               # Power input
        
""" Gas """   
        
class Gas:
    def __init__(self,size):
        self.size = size
        
    def set_spec(self,Esizemin,Esizemax, eps, lifetime):
        self.Esizemax=Esizemax              # [Wh] maximum storage size
        self.Esizemin=Esizemin              # [Wh] minimum storage size
        self.eps=eps                        # [MWh/MWh_stored] Energy cost per energy unit : 100 <=> it cost 100 J to
        self.lifetime = lifetime*8760
        
    def set_opti(self, capacity, energy):
        self.capacity = capacity            # Optimal storage capacity
        self.energy = energy                # Energy stored
        
