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

    print('sklearn : ' + str(pkg1))

    try:
        import networkx

        pkg5 = 'Pass'
    except:
        pkg5 = 'Fail'

    print('Networkx : ' + str(pkg5))