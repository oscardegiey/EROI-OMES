import pickle

def save_data(cells,grid,filename):

    # Saving the objects:
    with open('RESULTS/'+filename, 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([cells,grid], f)
        f.close()
    
def load_data(filename):
    
    # Getting back the objects:
    with open('RESULTS/'+filename,'rb') as f:  # Python 3: open(..., 'rb')
        [cells,grid] = pickle.load(f)
        f.close()
    
    return [cells,grid]
