from .dock import *
from .receptor import *
import time

class ScreenLibrary():
    
    
    def __init__(self,lig_set, name = 'Test'):
        
         #screen info
        self.name = name
        self.lib_path = str('./' + self.name + '/')
        
        #library to be screened
        self.dock_set = []
        self.rec_set = []
        self.lig_set = lig_set
        self.conf_set = []
        for i in range(len(lig_set)):
            self.conf_set += ([(conf , lig_set.loc[i, 'ligid']) for conf in lig_set.loc[i , 'ligand'].conformer_set])
        
        #files for input data
        self.rec_file = None       
   

    ##  ###  ## Determination of library contents ##  ###  ##
    
       
    #import receptor file funtion (called by collect_file())        
    def import_rec_file(self):
    
        rec_file = open(self.rec_file,'r')
        lines = rec_file.readlines()
        
        for line in lines:

            tabs = line.split()
            rec = Receptor()
        
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
        
    def create_dock_set(self):
                
            for rec in self.rec_set:
                
                c = 0
                for conf in self.conf_set:
                    
                    conf_id = str('l' + str(conf[1]) + 'c' + str(c))
                    
                    doc = Dock(rec, conf, conf_id)

                    self.dock_set.append(doc)
                    c += 1 
    
    
    ##  ###  ## Creation of library ##  ###  ##
    
    
    def prepare_library_directory_tree(self):
        
        try:
            #linux command (I dont think so anymore)
            try:
        
                sb.run('mkdir ' + self.name, shell = True)
            
                sb.run('mkdir ' + self.lib_path + 'pdb', shell = True)
                sb.run('mkdir ' + self.lib_path + 'pdbqt', shell = True)
                sb.run('mkdir ' + self.lib_path + 'res_pdbqt', shell = True)
                sb.run('mkdir ' + self.lib_path + 'res_sdf', shell = True)
                
                prep_pass = True
        
            #mac command??
            except:
            
                sb.run('mkdir ' + self.name, shell = True, stdout = sb.PIPE , text=True , check = True)
            
                sb.run('mkdir ' + self.lib_path + 'pdb', shell = True, stdout = sb.PIPE , text=True , check = True)
                sb.run('mkdir ' + self.lib_path + 'pdbqt', shell = True, stdout = sb.PIPE , text=True , check = True)
                sb.run('mkdir ' + self.lib_path + 'res_pdbqt', shell = True, stdout = sb.PIPE , text=True , check = True)
                sb.run('mkdir ' + self.lib_path + 'res_sdf', shell = True, stdout = sb.PIPE , text=True , check = True)
            
                prep_pass = True
            
        except:
            
            prep_pass = False
            print('\nDirectory Tree Creation Error\n')
                  
        
            
    def create_pdb_library(self):
        i = 0
        for conf in self.conf_set:
    
            AllChem.MolToPDBFile(conf[0], str(self.lib_path + 'pdb/l' + str(conf[1]) + 'c' + str(i) + '.pdb'))
            i += 1
        
    
    def create_pdbqt_library(self):
        
        problematic_ligands = 0
        
        i = 0
        for conf in self.conf_set:
            
            pdb_file = str(self.lib_path + 'pdb/l' + str(conf[1]) + 'c' + str(i) + '.pdb')
            pdbqt_file = str(self.lib_path + 'pdbqt/l' + str(conf[1]) + 'c' + str(i) + '.pdbqt')
                
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
                
            
            i += 1
            
        print('\nNumber of conversion errors: ' + str(problematic_ligands) + '\n')
            
    
    
    ##  ###  ## Screen ##  ###  ##
    
    
    
    def screen_library(self, ex = 8):
        
        problematic_ligands = 0
            
        for dock in self.dock_set:
            
            try:
                
                dock.vina_dock(self.lib_path, ex)
                
            except:
                
                problematic_ligands += 1
                
        print('\nNumber of docking errors: ' + str(problematic_ligands) + '\n')
        
    def convert_results(self):
        
        problematic_ligands = 0
        
        for dock in self.dock_set:
            
            in_file = str(self.lib_path + 'res_pdbqt/' + dock.name + '.pdbqt')
            out_file = str(self.lib_path + 'res_sdf/' + dock.name + '.sdf')
            
            command = str('obabel ' + in_file + ' -O' + out_file)
            
            try:
                
                sb.run(command, shell = True)
                
            except:
                
                problematic_ligands += 1
        
        print('\nNumber of conversion errors: ' + str(problematic_ligands) + '\n')
        
    def extract_results(self):
        
        problematic_ligands = 0 
        
        for dock in self.dock_set:
            
            file = str(self.lib_path + 'res_sdf/' + dock.name + '.sdf')
            
            mols = [mol for mol in AllChem.SDMolSupplier(file)]
            
            if mols is None:
                
                problematic_ligands += 1
            
            dock.conformations = mols
            
        print('\nNumber of empty results: ' + str(problematic_ligands) + '\n')
            
        ## NEXT ##
        #extract affinities
        
            #add (aff, mol) to dock.conformations
            #add 
   
    ##  ###  ## Clean Up ##  ###  ##
