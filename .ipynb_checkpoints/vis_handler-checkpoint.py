'''
module: vis_handler

purpose: functions for handling 3d visualizations in jupyter (possibly also image genreation depending
on output format)
'''

import py3Dmol
from IPython.display import display
from rdkit.Chem import AllChem
import ipywidgets

#ttps://birdlet.github.io/2019/10/02/py3dmol_example/

def show_mol3d(mol):
    
    #initiates a window of the specified size
    viewer = py3Dmol.view(600,300)
    
    viewer.addModel(AllChem.MolToMolBlock(mol) , 'mol') 
    
    viewer.setStyle({'stick' : {}})
    viewer.setViewStyle({'style':'outline','color':'black','width':0.1})
    
    #centres the molecule
    viewer.zoomTo()
    
    return viewer

def show_mols3d(mols):
    
    #initiates a window of the specified size
    viewer = py3Dmol.view(600,300)
    
    #add the molecule to the window
    for mol in mols:
        viewer.addModel(AllChem.MolToMolBlock(mol) , 'mol')
        
    #view options
    viewer.setStyle({'stick': {}})
    
    #centres the molecule
    viewer.zoomTo()
    
    return viewer

def show_dock_set(dock_set):
    
    def f(m, c): 
        return show_mol3d(dock_set[m].conformations[c])
    
    mol_slider = ipywidgets.IntSlider(min=0, max=len(dock_set)-1, step=1, description = 'Molecule')
    conf_slider = ipywidgets.IntSlider(min=0, max=8, step=1, description = 'Conformer') #max value needs fixing

    window = ipywidgets.interact(f,
                        m = mol_slider,
                        c = conf_slider)
    
    return window
    