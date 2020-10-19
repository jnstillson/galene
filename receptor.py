#!/usr/bin/env python
# coding: utf-8

class Receptor():
    
    def __init__(self):
        
        self.pos = [0.0, 0.0, 0.0]
        self.dim = [0.0, 0.0, 0.0]
        
        self.name = 'name' #must be unambiguous
        self.prot_name = 'protein'
        self.site_name = 'binding site'
        
        
        self.rig_file = 'rigid_file'
        self.flex_file = 'flex_file'
        
        self.flex = False





