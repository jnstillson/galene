#!/usr/bin/env python
# coding: utf-8

from .ligand import *
import time
import subprocess as sb


class Dock():

    def __init__(self, rec, conf, lig_id, conf_id):

        # parameters
        self.rec = rec
        self.conf = conf
        self.lig_id = lig_id
        self.conf_id = conf_id

        self.conf_tag = str('l' + str(lig_id) + 'c' + str(conf_id))
        self.name = str(rec.name + self.conf_tag)

        # debug
        self.wall_time = 0.0
        self.dock_pass = False

        self.vina_result_matrix = pd.DataFrame()
        self.vina_result_scores = pd.DataFrame()

        self.conformations = []

    def extract_conformations(self):

        # Maybe move conformer extraction here?

        # read in conformations as rdkit molecules and match them to their computed affinities
        pass

    def vina_dock(self,
                  lib_path,  # path where all files are stored
                  ex):  # exhaustiveness parameter

        vina_file = str(lib_path + 'res_pdbqt/' + self.name + '.pdbqt')
        conf_file = str(lib_path + 'pdbqt/' + self.conf_tag + '.pdbqt')
        time0 = time.time()

        dock_command = str('vina --out ' + vina_file +
                           ' --receptor ' + self.rec.rig_file +
                           ' --ligand ' + conf_file +
                           ' --center_x ' + str(self.rec.pos[0]) +
                           ' --center_y ' + str(self.rec.pos[1]) +
                           ' --center_z ' + str(self.rec.pos[2]) +
                           ' --size_x ' + str(self.rec.dim[0]) +
                           ' --size_y ' + str(self.rec.dim[1]) +
                           ' --size_z ' + str(self.rec.dim[2]) +
                           ' --exhaustiveness ' + str(ex))
        # print('\nDocking to receptor: '+str(self.rec.name)+'\n')

        try:

            sb.run(dock_command, shell=True)

            self.dock_pass = True

        except:

            self.dock_pass = False

            print('\nsomething happened ¯\_(ツ)_/¯\n')

        self.wall_time = time.time() - time0

        # print(str('this doc took ' + str(self.wall_time) + ' long.'))
