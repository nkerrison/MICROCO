# -------------------------------------------------------------------IMPORTS--------

from flask import Flask, render_template, request, send_file
import numpy as np
import pandas as pd
import re

# ----------------------------------------------------------------------------------

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

# ----------------------------------------------------------------------------------

@app.route('/upload', methods=['POST'])
def upload():
   
    # --- Request Excel file from the HTML form
    file   = request.files['file']
    oxprop = request.form.get('oxprop')
    oxprop = int(oxprop)

    # --- Read the Excel file and create separate dataframes of necessary sheets
    user_data_elox = pd.read_excel(file, sheet_name='El-Ox')
    user_data_stat = pd.read_excel(file, sheet_name='Stat')

    # --- Find the index (location) of the cell 'Weight%'
    weight_index = user_data_elox.columns.get_loc('Weight%')

    # --- Find the index (location) of the cell 'Oxide'
    oxide_index = user_data_elox.columns.get_loc('Oxide')

    # --- Extract columns pertaining to 'Weight%' (needed for Det Lim)
    weight_columns = user_data_elox.iloc[:, weight_index:oxide_index]

    # --- Extract columns pertaining to 'Oxide'
    oxide_columns = user_data_elox.iloc[:, oxide_index:]

    # --- Set the second row (the oxide names) as the header
    oxide_columns = oxide_columns.iloc[1:].rename(columns=oxide_columns.iloc[0])

    # --- Remove any rows that are entirely made of zeros to avoid
    # --- a division by zero error. Rows of zeros can occur, e.g.,
    # --- when acquision of a data point is canceled.
    oxide_columns_print = oxide_columns.loc[(oxide_columns.sum(axis=1)!=0) | (oxide_columns.applymap(lambda x: isinstance(x,str))).any(axis=1)]

    # --- Removes the last column ("Total") from the dataframe
    oxide_columns = oxide_columns_print.iloc[: , :-1]

    # ----------------------------------------------------------------------------------

    def get_formula_mass(formula):

        # --- Abridged atomic weights via Prohaska et al. 2022, "Standard atomic weights
        # --- of the elements 2021 (IUPAC Technical Report)", Table 1

        abridged_atomic_weights = {
              'H': 1.0080,   'He': 4.0026,   'Li': 6.94,     'Be': 9.0122,
              'B': 10.81,     'C': 12.011,    'N': 14.007,    'O': 15.999,
              'F': 18.998,   'Ne': 20.180,   'Na': 22.990,   'Mg': 24.305,
             'Al': 26.982,   'Si': 28.085,    'P': 30.974,    'S': 32.06,
             'Cl': 35.45,     'K': 39.098,   'Ar': 39.95,    'Ca': 40.078,
             'Sc': 44.956,   'Ti': 47.867,    'V': 50.942,   'Cr': 51.996,
             'Mn': 54.938,   'Fe': 55.845,   'Ni': 58.693,   'Co': 58.933,
             'Cu': 63.546,   'Zn': 65.38,    'Ga': 69.723,   'Ge': 72.630,
             'As': 74.922,   'Se': 78.971,   'Br': 79.904,   'Kr': 83.798,
             'Rb': 85.468,   'Sr': 87.62,     'Y': 88.906,   'Zr': 91.224,
             'Nb': 92.906,   'Mo': 95.95,    'Tc': 98.0,     'Ru': 101.07,
             'Rh': 102.91,   'Pd': 106.42,   'Ag': 107.87,   'Cd': 112.41,
             'In': 114.82,   'Sn': 118.71,   'Sb': 121.76,    'I': 126.90,
             'Te': 127.60,   'Xe': 131.29,   'Cs': 132.91,   'Ba': 137.33,
             'La': 138.91,   'Ce': 140.12,   'Pr': 140.91,   'Nd': 144.24,
             'Pm': 145.0,    'Sm': 150.36,   'Eu': 151.96,   'Gd': 157.25,
             'Tb': 158.93,   'Dy': 162.50,   'Ho': 164.93,   'Er': 167.26,
             'Tm': 168.93,   'Yb': 173.05,   'Lu': 174.97,   'Hf': 178.49,
             'Ta': 180.95,    'W': 183.84,   'Re': 186.21,   'Os': 190.23,
             'Ir': 192.22,   'Pt': 195.08,   'Au': 196.97,   'Hg': 200.59,
             'Tl': 204.38,   'Pb': 207.2,    'Bi': 208.98,   'Th': 232.04,
             'Pa': 231.04,    'U': 238.03,   'Np': 237.0,    'Pu': 244.0,
             'Am': 243.0,    'Cm': 247.0,    'Bk': 247.0,    'Cf': 251.0,
             'Es': 252.0,    'Fm': 257.0,    'Md': 258.0,    'No': 259.0,
             'Lr': 262.0,    'Rf': 267.0,    'Db': 270.0,    'Sg': 271.0,
             'Bh': 270.0,    'Hs': 277.0,    'Mt': 276.0,    'Ds': 281.0,
             'Rg': 280.0,    'Cn': 285.0,    'Nh': 284.0,    'Fl': 289.0,
             'Mc': 288.0,    'Lv': 293.0,    'Ts': 294.0,    'Og': 294.0
        }

        # --- Given a string, extract element symbols and the number (if any) that follows
        elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)

        formula_mass = 0
        
        for element, count in elements:
            
            atomic_weight = abridged_atomic_weights.get(element, 0)
            
            if count:
                
                formula_mass += int(count) * atomic_weight
            
            else:
                
                formula_mass += atomic_weight
        
        return formula_mass
    
    # ----------------------------------------------------------------------------------

    def count_oxygen_atoms(chemical_formula):

        # --- Given a string, return the number that follows the letter O
        oxygen_matches = re.findall(r'O(\d*)', chemical_formula)

        # --- Initialize the count
        oxygen_count = 0

        for match in oxygen_matches:

            # --- If a number succeeds the letter O, add to count
            if match:
                
                oxygen_count += int(match)

            # --- If no number succeeds the letter O, assume count is 1
            else:
                
                oxygen_count += 1

        return oxygen_count
    
    # ----------------------------------------------------------------------------------

    def ratio_of_cation_to_oxygen(oxide_formula):

        pattern = r'([A-Z][a-z]*)(\d*)'

        matches = re.findall(pattern, oxide_formula)

        cation_count = 0
        oxygen_count = 0

        for element, count in matches:

            if element == 'O':
                
                oxygen_count = int(count) if count else 1

            else:
                
                cation_count += int(count) if count else 1

        ratio = cation_count / oxygen_count

        return ratio
    
    # ----------------------------------------------------------------------------------

    def oxide_to_cation_str_converter(column_name):
        
        match = re.search(r'^[A-Z][a-z]*', column_name)
        
        if match:
            
            return match.group()
        
        else:
            
            return column_name

    # ----------------------------------------------------------------------------------

    # --- Grab the formula masses for each oxide in order to divide each
    # --- oxide weight by its corresponding formula mass
    formula_masses = [get_formula_mass(x) for x in oxide_columns.columns]
    molecular_proportions = oxide_columns.div(formula_masses, axis=1)

    # --- To derive the atomic proportion of oxygen from each molecule,
    # --- we must multiply the molecular proportion of oxide_columns by the
    # --- number of oxygen atoms in the oxide concerned.
    oxygen_counts = [count_oxygen_atoms(x) for x in oxide_columns.columns]
    atomic_proportions = molecular_proportions.mul(oxygen_counts, axis=1)

    # --- Recast oxygen atom proportions to the user's input
    tmp = pd.DataFrame()
    tmp['T2'] = oxprop / atomic_proportions.sum(axis=1)
    number_anions = atomic_proportions.mul(tmp['T2'], axis=0)

    # --- Calculate cation fractions
    ratios = [ratio_of_cation_to_oxygen(x) for x in oxide_columns.columns]
    cation_fractions = number_anions.mul(ratios, axis=1)
    cation_fractions['Total'] = cation_fractions.sum(axis=1)

    # --- Apply the function to all column names to extract prefixes
    prefixes = cation_fractions.columns.map(oxide_to_cation_str_converter)

    # --- Replace the column names with the extracted prefixes
    cation_fractions.columns = prefixes

    # --- Render template with DataFrames
    return render_template('result.html',
                           original_tables=[oxide_columns_print.to_html(classes='data')], 
                           processed_tables=[cation_fractions.to_html(classes='data')], 
                           titles=cation_fractions.columns.values)

if __name__ == '__main__':
    app.run(debug=True)