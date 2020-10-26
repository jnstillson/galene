from flask import Flask, render_template,redirect,url_for,send_file
from scripts.file_handler import *
from rdkit.Chem import Draw
from scripts.screen_library import ScreenLibrary
import io

#Variable Handling for Pages
lig_set = lig_set_from_smi('./data/lig/test_ligs.smi')
rec_set = rec_set_from_file('./data/rec2.txt')

screen = ScreenLibrary(lig_set, name = 'data/threeLigTest')

screen.rec_file = '/Users/jake/Desktop/galene/data/rec2.txt'
screen.import_rec_file()
screen.create_dock_set()
screen.extract_results()#why wont this display?
print(screen.lig_set.loc[0,'ligand'].vina_result_scores[0])

#Pages
app = Flask(__name__)

#MAIN PAGES
@app.route('/')
@app.route('/main')
def mainPage():
    return render_template('mainpage.html')

@app.route('/ligset')
def ligSetPage():
    return render_template('ligset.html', lig_set = lig_set)

@app.route('/recset')
def recSetPage():
    return render_template('recset.html', rec_set = rec_set)

@app.route('/progress')
def progressPage():
    return render_template('progress.html')

#OTHER DIRECTORIES
def pilImage(pil):
    img = io.BytesIO()
    pil.save(img, 'JPEG', quality = 300)
    img.seek(0)
    return send_file(img, mimetype='image/jpeg')

@app.route(f'/lig/<name>')
def ligPage(name):
    #look for lig name (should it be name?)
    if name in list(lig_set['name']):
        i = list(lig_set['name']).index(name)
        lig = lig_set.loc[i, 'ligand']
        img = Draw.MolToImage(lig.mol)
        return render_template('lig.html',
                               name=name,
                               lig=lig,
                               img=pilImage(img))
    else:
        return redirect(url_for('ligSetPage'))

@app.route(f'/lig/img/<name>.jpeg')
def ligImg(name):
    #look for lig name (should it be name?)
    if name in list(lig_set['name']):
        i = list(lig_set['name']).index(name)
        lig = lig_set.loc[i, 'ligand']
        img = Draw.MolToImage(lig.mol)
        return pilImage(img)
    else:
        return redirect(url_for('ligSetPage'))


@app.route(f'/rec/<name>')
def recPage(name):
    #for for rec name
    if name in list(rec_set['name']):
        i = list(rec_set['name']).index(name)
        rec = rec_set.loc[i, 'receptor']
        return render_template('rec.html', name=name, rec=rec)
    else:
        return redirect(url_for('recSetPage'))


if __name__ == '__main__':
    app.run(debug=True)

