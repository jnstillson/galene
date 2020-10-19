#!/usr/bin/env python
# coding: utf-8

from .ligand import *
from .receptor import *
import time
import IPython
import subprocess as sb
import sys
import os

class Dock():
    
    def __init__(self, rec, conf, conf_id):
        
        #parameters
        self.rec = rec
        self.conf = conf
        self.conf_id = conf_id
        
        self.name = str(rec.name[:4] + conf_id)

        #debug
        self.wall_time = 0.0
        self.dock_pass = False
        
        #results
        self.conformations = []
        
    def extract_conformations(self):
        
        # Maybe move conformer extraction here?
        
        # read in conformations as rdkit molecules and match them to their computed affinities
        pass    
    
    def vina_dock(self,
                  lib_path, #path where all files are stored
                  ex): # exhaustiveness parameter
        
        vina_file = str(lib_path + 'res_pdbqt/'+ self.name + '.pdbqt')
        conf_file = str(lib_path + 'pdbqt/' + self.conf_id + '.pdbqt')
        
        time0 = time.time()
        
    
        dock_command = str('vina --out '+ vina_file + 
                   ' --receptor ' + self.rec.rig_file +
                   ' --ligand ' + conf_file + 
                   ' --center_x ' + str(self.rec.pos[0]) + 
                   ' --center_y ' + str(self.rec.pos[1]) + 
                   ' --center_z ' + str(self.rec.pos[2]) +
                   ' --size_x ' + str(self.rec.dim[0]) + 
                   ' --size_y ' + str(self.rec.dim[1]) +
                   ' --size_z ' + str(self.rec.dim[2]) +
                   ' --exhaustiveness ' + str(ex))
                
        
        if self.rec.flex == True:
            
            dock_command += str(' --flex ' + self.rec.flex_file)
    
        #print('\nDocking to receptor: '+str(self.rec.name)+'\n')
    
        try:
            
            sb.run(dock_command, shell = True)
        
            self.dock_pass = True
            
        except:
            
            self.dock_pass = False
            
            print('\nsomething happened ¯\_(ツ)_/¯\n')
            
        
        self.wall_time = time.time() - time0
        
        #print(str('this doc took ' + str(self.wall_time) + ' long.'))
