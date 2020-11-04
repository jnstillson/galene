from flask import Flask, render_template,redirect,url_for,send_file, request
from scripts.file_handler import *
from scripts.util import *
from rdkit.Chem import Draw
from scripts.screen_library import ScreenLibrary
from scripts.analyze_lig_set import AnalyzeLigSet
from scripts.analyze_dock_poses import *
import io
import subprocess as sb

global lig_set
global rec_set
global display_set

display_set = {'lig_set': lig_set_from_smi('./data/cofactors.smi'),
               'rec_set': rec_set_from_file('./data/rec2.txt'),
               'pca_plot': None,
               'frag_plot': None,
               'frag_set': None,
               }

display_set['pca_plot'] = AnalyzeLigSet(lig_set = display_set['lig_set'], pca_plot = 1).pca_plot
display_set['lig_set_prop'] = get_lig_set_properties(display_set['lig_set'])
display_set['rec_set_prop'] = get_rec_set_properties(display_set['rec_set'])
#Variable Handling for Pages

lig_set = display_set['lig_set']

#Pages
app = Flask(__name__)

#MAIN PAGES
@app.route('/')
@app.route('/main')
@app.route('/home')
def mainPage():
    return render_template('home.html')

@app.route('/ligset')
def ligSetPage():
    display_set['lig_set_prop'] = get_lig_set_properties(display_set['lig_set'])
    return render_template('ligset.html', lig_set = display_set['lig_set'], props = display_set['lig_set_prop'])

@app.route('/recset')
def recSetPage():
    display_set['rec_set_prop'] = get_rec_set_properties(display_set['rec_set'])
    return render_template('recset.html', rec_set = display_set['rec_set'], props = display_set['rec_set_prop'])

@app.route('/progress')
def progressPage():
    return render_template('progress.html')

#OTHER DIRECTORIES
def pilImage(pil):
    img = io.BytesIO()
    pil.save(img, 'JPEG', quality = 300)
    img.seek(0)
    return send_file(img, mimetype='image/jpeg')

@app.route('/rec/data/rec/pdb/<pdb>')
def pdbFiles(pdb):
    for rec in display_set['rec_set']['receptor']:
        return send_file(rec.pdb_file)


@app.route(f'/lig/<name>', methods= ['GET', 'POST'])
def ligPage(name):
    recid = 0
    if request.method == 'POST':
        for i in range(len(display_set['rec_set']['receptor'])):
            if request.form['rec'] == display_set['rec_set']['name'][i]:
                recid = display_set['rec_set']['recid'][i]
    lig_set = display_set['lig_set']
    if name in list(lig_set['name']):
        i = list(lig_set['name']).index(name)
        lig = lig_set.loc[i, 'ligand']
        img = Draw.MolToImage(lig.mol)
        return render_template('lig.html',
                               name=name,
                               lig=lig,
                               img=pilImage(img),
                               molBlock=lig.mol_block,
                               display_set = display_set,
                               recid = recid,
                               )
    else:
        return redirect(url_for('ligSetPage'))

@app.route(f'/lig/img/<name>.jpeg')
def ligImg(name):
    lig_set = display_set['lig_set']
    if name in list(lig_set['name']):
        i = list(lig_set['name']).index(name)
        lig = lig_set.loc[i, 'ligand']
        img = Draw.MolToImage(lig.mol)
        return pilImage(img)
    else:
        return redirect(url_for('ligSetPage'))


@app.route(f'/rec/<name>')
def recPage(name):
    rec_set = display_set['rec_set']
    if name in list(rec_set['name']):
        i = list(rec_set['name']).index(name)
        rec = rec_set.loc[i, 'receptor']
        return render_template('rec.html', name=name, rec=rec)
    else:
        return redirect(url_for('recSetPage'))

@app.route('/setup', methods=['POST','GET'])
def setupPage():
    ligSuc = 'No File'
    recSuc = 'No File'
    resSuc = 'No File'
    if request.method == 'POST':
        if request.form['upload'] == 'Upload':
            ligFile = request.values.get('lig file')
            recFile = request.values.get('rec file')
            resFolder = request.values.get('result folder')
            print(ligFile)
            print(recFile)
            print(resFolder)
            try:
                display_set['rec_set'] = rec_set_from_file(recFile)
                recSuc = recFile
            except:
                recSuc = 'No New File'
            if ligFile[-3:] == 'smi':
                print('hey')
                try:
                    display_set['lig_set'] = lig_set_from_smi(ligFile)
                    ligSuc = ligFile
                except:
                    ligSuc = 'No New File'
            elif ligFile[-3:] =='sdf':
                try:
                    display_set['lig_set'] = lig_set_from_gypsum(ligFile)
                    ligSuc = ligFile
                except:
                    ligSuc = 'No New File'
            else:
                ligSuc = 'No New File'


            try:
                screen = ScreenLibrary(display_set['lig_set'], display_set['rec_set'], name=resFolder)
                screen.create_dock_set()
                screen.extract_results()
                resSuc = resFolder
            except:
                resSuc = 'No New Results'

    return render_template('setup.html', ligSuc=ligSuc, recSuc=recSuc, resSuc=resSuc)

@app.route('/dock', methods = ['GET', 'POST'])
def dockPage():
    display_set['lig_set_prop'] = get_lig_set_properties(display_set['lig_set'])
    display_set['rec_set_prop'] = get_rec_set_properties(display_set['rec_set'])
    ex = 1
    name = 'Name'
    if request.method == 'POST':
        if request.form['dock'] == 'Update Parameters':
            ex = int(request.values.get('ex'))
            name = str(request.values.get('name'))
            #update commands go here
        if request.form['dock'] == 'Start Dock':
            ex = int(request.values.get('ex'))
            name = str(request.values.get('name'))
            if check_name(name):
                try:
                    screen = ScreenLibrary(display_set['lig_set'],
                                       display_set['rec_set'],
                                       name= name,)
                    screen.create_dock_set()
                    screen.prepare_library_directory_tree()
                    screen.create_pdb_library()
                    screen.create_pdbqt_library()
                    screen.screen_library(ex = ex)
                    screen.convert_results()
                    screen.extract_results()
                except:
                    pass
            else:
                name = 'Bad Name'
            #dock commands go here
    return render_template('dock.html', display_set = display_set, ex = ex, name=name)

@app.route('/analyze', methods = ['GET', 'POST'])
def analyzePage():
    cluster = 0
    fcluster = 0
    if request.method == 'POST':
        if request.form['recalc'] == 'Recalculate Ligand Distribution':
            cluster = request.form['clustNum']

            pca = AnalyzeLigSet(display_set['lig_set'], pca_plot=cluster)
            display_set['pca_plot'] = pca.pca_plot
            display_set['pca_raw'] = pca.lig_space2
            #call analyzeligset with new pca_plot method that takes in an int.
        if request.form['recalc'] == 'Recalculate Fragment Distribution':
            fcluster = request.form['fclustNum']
            display_set['frag_plot'] = AnalyzeLigSet(display_set['lig_set'], frag_plot=fcluster).fpca_plot

    return render_template('analyze.html',
                           ligPlot = display_set['pca_plot'],
                           fragPlot = display_set['frag_plot'],
                           cluster = cluster,
                           fcluster = fcluster,
                           )

@app.route('/results', methods=['GET', 'POST'])
def resultsPage():

    if request.method == 'POST':
        rec = request.form['rec']
        id = display_set['rec_set'][display_set['rec_set']['name'] == rec]['recid'].values[0]
        print(id)
        best_scores = min_score(display_set['lig_set'], id = id)
        colors = calculate_colors(best_scores)
        display_set['result_plot'] = scatter_to_js_color(display_set['pca_raw'], colors)

    return render_template('results.html', display_set = display_set)

@app.route('/settings')
def settingsPage():
    return render_template('settings.html')


if __name__ == '__main__':
    app.run(debug=True)

