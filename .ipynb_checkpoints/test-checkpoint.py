#### ## #### Testing #### ## ####

## create and prep smaller screen
test_small = screen_library('test_small')

test_small.rec_file = '/Users/jake/Desktop/paid/rec2.txt'

test_small.lig_file = '/Users/jake/Desktop/paid/cofactors.smi'

test_small.import_rec_file()

test_small.import_lig_file()

test_small.create_dock_set()

test_small.prepare_library_directory_tree()

test_small.create_pdb_library()

test_small.create_pdbqt_library()

test_small.screen_library()

test_small.convert_results()