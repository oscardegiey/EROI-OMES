
class Grid:
    """Calss representing the grid.
    
        Attributes:
            e_imp_max: maximum electricity import in a given cell
            e_exch_max: maximum electricity exchange between 2 given cells
            gas_imp_max: maximum gas import in a given cell
            gas_exch_max: maximum gas exchange between 2 given cells
            
            e_import: vector of imported electricity in a cell at each time t
            e_exch: vector of exchanged electricity between 2 cells at each time t
            gas_import: vector of imported electricity in a cell at each time t
            gas_exch: vector of exchanged gas between 2 cells at each time t 
            
    """
    
    def __init__(self, e_imp_max, e_exch_max, gas_imp_max, gas_exch_max):
        """Inits Grid."""
        self.e_imp_max = e_imp_max
        self.e_exch_max = e_exch_max
        self.gas_imp_max = gas_imp_max
        self.gas_exch_max = gas_exch_max
        
    def set_opti(self, e_import, e_exch, gas_import, gas_exch):
        """Set the optimal value."""
        self.e_import = e_import
        self.e_exch = e_exch
        self.gas_import = gas_import
        self.gas_exch = gas_exch

    