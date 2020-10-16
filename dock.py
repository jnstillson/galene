#!/usr/bin/env python
# coding: utf-8

from .ligand import *
from .receptor import *
import time
import IPython
import subprocess as sb
import sys
import os

class dock():
    
    def __init__(self):
        
        #parameters
        self.rec = receptor()
        self.lig = ligand()
        
        self.name = str(self.rec.name[:4] + '_' + self.lig.name[:4])

        #debug
        self.wall_time = 0.0
        self.dock_pass = False
        
        #results
        self.conformations = []
        
    def extract_conformations(self):
        
        # read in conformations as rdkit molecules and match them to their computed affinities
        pass    
    
    def vina_dock(self, 
                  vina_file, # path for .pdbqt output file
                  lig_file,  # path for .pdbqt ligand input file
                  rec_file,  # path for .pdbqt (rigid) receptor input file
                  flex_file, # path for .pdbqt flexible receptor input file
                  exhaustivness = 1 # exhaustiveness parameter
                 ): 
    
        time0 = time.time()
    
        dock_command = str('vina --out '+ vina_file + 
                   ' --receptor ' + rec_file +
                   ' --ligand ' + lig_file + 
                   ' --center_x ' + str(self.rec.pos[0]) + 
                   ' --center_y ' + str(self.rec.pos[1]) + 
                   ' --center_z ' + str(self.rec.pos[2]) +
                   ' --size_x ' + str(self.rec.dim[0]) + 
                   ' --size_y ' + str(self.rec.dim[1]) +
                   ' --size_z ' + str(self.rec.dim[2]) +
                   ' --exhaustiveness ' + str(exhaustivness))
                
        if self.rec.flex == True:
            
            dock_command += str(' --flex ' + flex_file)
    
        print('\nReceptor: '+str(self.rec.name)+'\n\nLigand: '+str(self.lig.name)+'\n')
    
        try:
            
            sb.run(dock_command, shell = True)
        
            self.dock_pass = True
            
        except:
            
            self.dock_pass = False
            
            print('\nsomething happened ¯\_(ツ)_/¯\n')
            
        
        self.wall_time = time.time() - time0
        
        print(str('this doc took ' + str(self.wall_time) + ' long.'))





