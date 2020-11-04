'''
This file will be for tracking utility stuff, wall times, error handling etc


'''
import numpy as np

def scatter_to_js(scatter, clusters):
    # assuming the structure of the plot is [[x vals],[y vals]]
    cluster_colors = ['#CFFBE2', '#AFDAF7', '#EFE583', '#B5AE4F',
                      '#F9CD9B', '#D3CDFF', '#E1D776', '#F2A694', ]
    scatter = np.array(scatter)
    scatter_js = '{ datasets: [ '
    for c in range(max(clusters) + 1):
        scatter_js += '{'
        scatter_js += f' label: "Cluster {c}", data: ['
        for i in range(len(scatter)):
            if clusters[i] == c:
                scatter_js += str('{ x: '
                                  + str(scatter[i][0])
                                  + ', y: '
                                  + str(scatter[i][1])
                                  + '},\n')
        scatter_js += f'], backgroundColor : "{cluster_colors[c]}"'
        scatter_js += '},'
    scatter_js += '],}'
    return scatter_js

def calculate_colors(scores, min_color="#A9C9FF", max_color="#F2A694"):
    colors = []
    max_score = -7.0
    min_score = -10.2
    for scr in scores:
        if scr[0] <= max_score:
            val = int(abs((scr[0]-max_score)/(max_score-min_score)) * 255)
        else:
            val = 0
        colors += [f'rgb({255-val}, {val}, 10)']

    return colors


def scatter_to_js_color(scatter, colors):
    scatter = np.array(scatter)
    scatter_js = '{ datasets: [{ label:"results", data: ['
    for i in range(len(scatter)):
        scatter_js += str('{ x: '
                            + str(scatter[i][0])
                            + ', y: '
                            + str(scatter[i][1])
                            + '},\n')
    scatter_js += '], backgroundColor : ['
    for i in range(len(colors)):
        scatter_js += f'"{colors[i]}",'

    scatter_js += ']}]}'
    return scatter_js


def get_lig_set_properties(lig_set):

    mmass = []
    logp = []
    mr = []
    h_don = []
    h_acc = []

    for lig in lig_set['ligand']:
        mmass.append(lig.mmass)
        logp.append(lig.logp)
        mr.append(lig.mr)
        h_don.append(lig.h_don)
        h_acc.append(lig.h_acc)

    lig_set_prop = {'size' : len(lig_set['ligand']),
                    'av mmass' : sum(mmass)/len(mmass),
                    'av logp' : sum(logp)/len(logp),
                    'av mr' : sum(mr)/len(mr),
                    'av h_don' : sum(h_don)/len(h_don),
                    'av h_acc' : sum(h_acc)/len(h_acc)
                    }
    return lig_set_prop

def get_rec_set_properties(rec_set):

    rec_set_prop = {'size' : len(rec_set),}

    return rec_set_prop

def check_name(name):

    #make sure its not empty
    if name is None:
        return False

    if name == '':
        return False

    #check for bad characters
    bad_char = [' ','.','/','?','*','\\','<','>','|',':','\"','\'']
    for b in bad_char:
        if b in name:
            return False

    #check to see if name is already being used

    return True


class Error(Exception):
    pass


# need screen_library debug


# analyze_lig_set debug
class LigSetAnalysisError(Error):

    def __init__(self, message='\nError while precomputing values for unsupervised ligand analysis\n'):
        self.message = message
