#### ## #### Testing #### ## ####
if __name__ == '__main__':
    try:
        import rdkit.Chem
        pkg1 = 'Pass'
    except:
        pkg1 = 'Fail'

    print('RDKit : ' + str(pkg1))

    try:
        import flask
        pkg2 = 'Pass'
    except:
        pkg2 = 'Fail'

    print('Flask : ' + str(pkg2))

    try:
        import pandas
        pkg3 = 'Pass'
    except:
        pkg3 = 'Fail'

    print('Pandas : ' + str(pkg3))

    try:
        import sklearn
        pkg4 = 'Pass'
    except:
        pkg4 = 'Fail'

    print('sklearn : ' + str(pkg4))

    try:
        import networkx
        pkg5 = 'Pass'
    except:
        pkg5 = 'Fail'

    print('Networkx : ' + str(pkg5))

    try:
        import ipywidgets
        pkg6 = 'Pass'
    except:
        pkg6 = 'Fail'

    print('ipywidgets : ' + str(pkg6))

    try:
        import py3Dmol
        pkg7 = 'Pass'
    except:
        pkg7 = 'Fail'

    print('py3Dmol : ' + str(pkg7))

    try:
        import numpy
        pkg8 = 'Pass'
    except:
        pkg8 = 'Fail'

    print('Numpy : ' + str(pkg8))
