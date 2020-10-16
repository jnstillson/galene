from .dock import *
from .receptor import *
import time

class screen_library():
    
    
    def __init__(self, name = 'library_screen'):
        
        #library to be screened
        self.dock_set = []
        self.rec_set = []
        self.lig_set = []
        
        #screen info
        self.name = name
        self.lib_path = str('./' + self.name + '/')
        
        #files for input data
        self.rec_file = 'receptor file'
        self.lig_file = 'ligand file'
        
   

    ##  ###  ## Determination of library contents ##  ###  ##
    
    
    
    #import receptor file funtion (called by collect_file())        
    def import_rec_file(self):
    
        rec_file = open(self.rec_file,'r')
        lines = rec_file.readlines()
        
        for line in lines:

            tabs = line.split()
            rec = receptor()
        
            rec.name = str(tabs[0])
            rec.pos[0] = float(tabs[1])
            rec.pos[1] = float(tabs[2])
            rec.pos[2] = float(tabs[3])
            rec.dim[0] = float(tabs[4])
            rec.dim[1] = float(tabs[5])
            rec.dim[2] = float(tabs[6])
            rec.rig_file = str(tabs[7])
            rec.flex_file = str(tabs[8])
            
            if str(tabs[8]) != 'null':
                rec.flex = True
            
            self.rec_set.append(rec)
        
        rec_file.close()
    
    
    ## Formatted to read .smi files 'smiles \t name' ##
    
    #import ligand file function (called by collect_file())
    def import_lig_file(self):
        
        problematic_ligands = 0

        lig_file = open(self.lig_file,'r')
        lines = lig_file.readlines()
    
        for line in lines:
            tabs = line.split()
        
            lig = ligand()
        
            lig.smiles = str(tabs[0])
            lig.name = str(tabs[1])
            
            ## structure prep
            try:
                
                lig.prep_struct()
            
            except:
                
                lig.flag_issue()
                problematic_ligands += 1
        
            self.lig_set.append(lig)
            
        lig_file.close()
        
        print('\nNumber of ligands in set: ' + str(len(self.lig_set))+ '\n')
        
        print('\nNumber of prep errors: ' + str(problematic_ligands) + '\n')
        
    def create_dock_set(self):

        for rec in self.rec_set:
            for lig in self.lig_set:
                
                doc = dock()
                
                doc.rec = rec
                doc.lig = lig
                
                doc.name = str(rec.name+'_'+lig.name)
                doc.file = str(self.name + '/' + doc.name + '.pdbqt')
            
                self.dock_set.append(doc)
    
    
    
    ##  ###  ## Creation of library ##  ###  ##
    
    
    
    def prepare_library_directory_tree(self):
        
        try:
            #linux command
            try:
        
                sb.run('mkdir ' + self.name, shell = True)
            
                sb.run('mkdir ' + self.lib_path + 'pdb', shell = True)
                sb.run('mkdir ' + self.lib_path + 'pdbqt', shell = True)
                sb.run('mkdir ' + self.lib_path + 'results', shell = True)
            
                prep_pass = True
        
            #mac command??
            except:
            
                sb.run('mkdir ' + self.name, shell = True, stdout = sb.PIPE , text=True , check = True)
            
                sb.run('mkdir ' + self.lib_path + 'pdb', shell = True, stdout = sb.PIPE , text=True , check = True)
                sb.run('mkdir ' + self.lib_path + 'pdbqt', shell = True, stdout = sb.PIPE , text=True , check = True)
                sb.run('mkdir ' + self.lib_path + 'results', shell = True, stdout = sb.PIPE , text=True , check = True)
            
                prep_pass = True
            
        except:
            
            prep_pass = False
            print('\nDirectory Tree Creation Error\n')
                  
        
            
    def create_pdb_library(self):
        
        for dock in self.dock_set:
    
            dock.lig.file_writer(str(self.lib_path + 'pdb/'))
        
            
    def create_pdbqt_library(self):
        
        problematic_ligands = 0
        
        for dock in self.dock_set:
            
            pdb_file = str(self.lib_path + 'pdb/' + dock.lig.name + '.pdb')
            pdbqt_file = str(self.lib_path + 'pdbqt/' + dock.lig.name + '.pdbqt')
                
            command = str('obabel ' + pdb_file + ' -O' + pdbqt_file + ' --partialcharge gasteiger')
            
            try:
                
                try:
                
                    #linux command
                    sb.run(command, shell = True)
                
                except:
                
                    #mac command??
                    sb.run(command, shell = True, stdout = sb.PIPE , text=True , check = True)
            
            except:
                
                dock.lig.flag_issue()
                problematic_ligands += 1
                
                
        print('\nNumber of conversion errors: ' + str(problematic_ligands) + '\n')
            
    
    
    ##  ###  ## Screen ##  ###  ##
    
    
    
    def screen_library(self):
        
        problematic_ligands = 0
            
        for dock in self.dock_set:
            
            vina_file = str(self.lib_path + 'results/'+ dock.name + '.pdbqt')
            lig_file = str(self.lib_path + 'pdbqt/' + dock.lig.name + '.pdbqt')
                
            try:
                
                dock.vina_dock(vina_file,
                                lig_file,
                                dock.rec.rig_file,
                                dock.rec.flex_file)
                    
                
            except:
                
                problematic_ligands += 1
                
        print('\nNumber of docking errors: ' + str(problematic_ligands) + '\n')
        
    def convert_results(self):
        
        problematic_ligands = 0
        
        for dock in self.dock_set:
            
            in_file = str(self.lib_path + 'results/' + dock.name + '.pdbqt')
            out_file = str(self.lib_path + 'results/' + dock.name + '.mol2')
            
            command = str('obabel ' + in_file + ' -O' + out_file)
            
            try:
                
                sb.run(command, shell = True)
                
            except:
                
                problematic_ligands += 1
        
        print('\nNumber of conversion errors: ' + str(problematic_ligands) + '\n')
        
    def extract_results(self):
        
        problematic_ligands = 0 
        
        for dock in self.dock_set:
            
            file = str(self.lib_path + 'results/' + dock.name + '.mol2')
            
            mol = AllChem.MolFromMol2File(file)
            
            if mol is None:
                
                problematic_ligands += 1
            
            dock.conformations.append(mol)
            
        print('\nNumber of empty results: ' + str(problematic_ligands) + '\n')
            
        ## NEXT ##
        #extract affinities
        
            #add (aff, mol) to dock.conformations
            #add 
   
    ##  ###  ## Clean Up ##  ###  ##
