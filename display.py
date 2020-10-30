from flask import Flask, render_template,redirect,url_for,send_file, request
from scripts.file_handler import *
from scripts.util import *
from rdkit.Chem import Draw
from scripts.screen_library import ScreenLibrary
from scripts.analyze_lig_set import AnalyzeLigSet
import io
import subprocess as sb

global lig_set
global rec_set
global display_set

display_set = {'lig_set': lig_set_from_smi('./data/lig/dnmt1_zinc.smi'),
               'rec_set': rec_set_from_file('./data/rec2.txt'),
               'pca_plot': None,
               'frag_plot': None,
               'frag_set': None,
               }

display_set['pca_plot'] = AnalyzeLigSet(lig_set = display_set['lig_set'], pca_plot = True).pca_plot
display_set['lig_set_prop'] = get_lig_set_properties(display_set['lig_set'])
display_set['rec_set_prop'] = get_rec_set_properties(display_set['rec_set'])
#Variable Handling for Pages

lig_set = display_set['lig_set']

screen = ScreenLibrary(lig_set, name = 'threeLigTest')

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

@app.route(f'/lig/<name>', methods= ['GET', 'POST'])
def ligPage(name):
    curRec = display_set['rec_set']['receptor'][0].name
    if request.method == 'POST':
        for rec in display_set['rec_set']['receptor']:
            if request.form['rec'] == rec.name:
                curRec = str(rec.name)
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
                               curRec = curRec,
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

    ligFile = request.values.get('ligFile')
    recFile = request.values.get('recFile')
    try:
        try:
            display_set['lig_set'] = lig_set_from_smi(ligFile)
            ligSuc = 'New File Uploaded!'
        except:
            ligSuc = 'No New File'
        try:
            display_set['rec_set'] = rec_set_from_file(recFile)
            recSuc = 'New File Uploaded!'
        except:
            recSuc = 'No New File'
    except:
        print('Error')
    return render_template('setup.html', ligSuc=ligSuc, recSuc=recSuc)

@app.route('/dock', methods = ['GET', 'POST'])
def dockPage():
    if request.method == 'POST':
        if request.form['dock'] == 'Update Parameters':
            display_set['lig_set_prop'] = get_lig_set_properties(display_set['lig_set'])
            display_set['rec_set_prop'] = get_rec_set_properties(display_set['rec_set'])
            print('update here')
            #update commands go here
        if request.form['dock'] == 'Start Dock':
            print('dock here')
            #dock commands go here
    return render_template('dock.html', display_set = display_set)

@app.route('/analyze', methods = ['GET', 'POST'])
def analyzePage():
    cluster = 0
    fcluster = 0
    if request.method == 'POST':
        if request.form['recalc'] == 'Recalculate Ligand Distribution':
            cluster = request.form['clustNum']
            display_set['pca_plot'] = AnalyzeLigSet(display_set['lig_set'], pca_plot=cluster).pca_plot

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

@app.route('/results')
def resultsPage():
    return render_template('results.html')

@app.route('/settings')
def settingsPage():
    return render_template('settings.html')


if __name__ == '__main__':
    app.run(debug=True)

