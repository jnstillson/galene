from .dock import *
from .receptor import *
import time

class ScreenLibrary():

    def __init__(self,
                 lig_set,
                 rec_set,
                 name='Test',):

        # screen info
        self.name = name
        self.lib_path = str('./data/' + self.name + '/')

        # library to be screened
        self.dock_set = []
        self.rec_set = rec_set
        self.lig_set = lig_set
        self.conf_set = []
        self.times = []
        for i in range(len(lig_set)):
            self.conf_set += ([(conf,
                                lig_set.loc[i, 'ligand'].conformer_set.index(conf),  # index of the conformer
                                lig_set.loc[i, 'ligid'],  # index of the ligand
                                ) for conf in lig_set.loc[i, 'ligand'].conformer_set])

        # files for input data
        self.rec_file = None

        ##  ###  ## Determination of library contents ##  ###  ##

    def create_dock_set(self):

        for rec in self.rec_set['receptor']:

            c = 0
            for conf in self.conf_set:
                doc = Dock(rec,
                           conf,
                           conf_id=conf[1],
                           lig_id=self.lig_set.loc[conf[2], 'ligid'],
                           )

                self.dock_set.append(doc)
                c += 1

                ##  ###  ## Creation of library ##  ###  ##

    def prepare_library_directory_tree(self):

        try:
            # linux command (I dont think so anymore)
            try:

                sb.run('mkdir ' + self.lib_path, shell=True)

                sb.run('mkdir ' + self.lib_path + 'pdb', shell=True)
                sb.run('mkdir ' + self.lib_path + 'pdbqt', shell=True)
                sb.run('mkdir ' + self.lib_path + 'res_pdbqt', shell=True)
                sb.run('mkdir ' + self.lib_path + 'res_sdf', shell=True)

                prep_pass = True

            # mac command??
            except:

                sb.run('mkdir ' + self.lib_path, shell=True, stdout=sb.PIPE, text=True, check=True)
                sb.run('mkdir ' + self.lib_path + 'pdb', shell=True, stdout=sb.PIPE, text=True, check=True)
                sb.run('mkdir ' + self.lib_path + 'pdbqt', shell=True, stdout=sb.PIPE, text=True, check=True)
                sb.run('mkdir ' + self.lib_path + 'res_pdbqt', shell=True, stdout=sb.PIPE, text=True, check=True)
                sb.run('mkdir ' + self.lib_path + 'res_sdf', shell=True, stdout=sb.PIPE, text=True, check=True)

                prep_pass = True

        except:

            prep_pass = False
            print('\nDirectory Tree Creation Error\n')

    def create_pdb_library(self):
        i = 0
        for conf in self.conf_set:
            AllChem.MolToPDBFile(conf[0], str(self.lib_path + 'pdb/l' + str(conf[2]) + 'c' + str(conf[1]) + '.pdb'))
            i += 1

    def create_pdbqt_library(self):

        problematic_ligands = 0

        i = 0
        for conf in self.conf_set:

            pdb_file = str(self.lib_path + 'pdb/l' + str(conf[2]) + 'c' + str(conf[1]) + '.pdb')
            pdbqt_file = str(self.lib_path + 'pdbqt/l' + str(conf[2]) + 'c' + str(conf[1]) + '.pdbqt')

            command = str('obabel ' + pdb_file + ' -O' + pdbqt_file + ' --partialcharge gasteiger')

            try:

                try:

                    # linux command
                    sb.run(command, shell=True)

                except:

                    # mac command??
                    sb.run(command, shell=True, stdout=sb.PIPE, text=True, check=True)

            except:

                problematic_ligands += 1

            i += 1

        print('\nNumber of conversion errors: ' + str(problematic_ligands) + '\n')

    ##  ###  ## Screen ##  ###  ##

    def screen_library(self, ex=8):

        problematic_ligands = 0

        for dock in self.dock_set:
            t0 = time.time()


            try:

                dock.vina_dock(self.lib_path, ex)
                docking_time = time.time()-t0
                self.times.append(docking_time)

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

                sb.run(command, shell=True)

            except:

                problematic_ligands += 1

        print('\nNumber of conversion errors: ' + str(problematic_ligands) + '\n')

    def extract_results(self, path = None):
        if path == None:
            path = self.lib_path
        problematic_ligands = 0
        for lig in self.lig_set['ligand']:
            lig.vina_result_scores = [pd.DataFrame() for _ in range(len(self.rec_set['receptor']))]
            lig.vina_result_poses = [pd.DataFrame() for _ in range(len(self.rec_set['receptor']))]

        for dock in self.dock_set:
            file = str(path + 'res_sdf/' + dock.name + '.sdf')
            mols = [mol for mol in AllChem.SDMolSupplier(file)]
            dock.conformations = mols

            # attach results to ligand object
            lig = self.lig_set.loc[dock.lig_id, 'ligand']
            lig.cheminformatic_analysis()
            r = list(self.rec_set['receptor']).index(dock.rec)
            for i in range(len(mols)):
                try:
                    score = mols[i].GetProp('REMARK').split()[2]
                    lig.vina_result_scores[r].loc[i, dock.conf_id] = float(score)
                    lig.vina_result_poses[r].loc[i, dock.conf_id] = mols[i]
                except:
                    problematic_ligands += 1

        print('\nNumber of empty results: ' + str(problematic_ligands) + '\n')

        #potentially change to be a single matrix per receptor

    def write_summary(self, description = None, ex = 0, ):
        summary = open(self.lib_path + 'screen_report.txt', 'w+')

        summary.write('### Screening Job Info ###\n\n')

        summary.write('Name:\t' + self.name + '\n')
        summary.write('Exhaustiveness:\t' + str(ex) + '\n')
        summary.write('Receptor Number:\t' + str(len(self.rec_set)) + '\n')
        summary.write('Ligand Number:\t' + str(len(self.lig_set)) + '\n')
        summary.write('Description:\t' + description + '\n')

        summary.write('\n### Screening Results ###\n\n')
        summary.write('Average Docking Time:\t' + str(sum(self.times)/len(self.times)) + '\n')

        summary.write('\n### Raw Results / Debug ###\n\n')
        summary.write(str(self.times))

        summary.close()
        print('done.')


    ##  ###  ## Clean Up ##  ###  ##
