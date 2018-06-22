import logging
import os
from pastml import F81, MARGINAL_APPROXIMATION

import pandas as pd
from ete3 import Tree
from hdx.location.country import Country

from cytopast import remove_certain_leaves
from cytopast.pastml_analyser import pastml_pipeline
from generate_geomap import generate_map

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, '4', 'pastml_fast_tree.nwk')
SA_SUBTREE_NWK = os.path.join(DATA_DIR, '4', 'pastml_fast_tree.sa.nwk')
YEAR_SUBTREE_NWK = os.path.join(DATA_DIR, '4', 'pastml_fast_tree.{}.nwk')
STATES_INPUT = os.path.join(DATA_DIR, 'metadata.tab')
STATES_INPUT_LOC = os.path.join(DATA_DIR, 'metadata_loc.tab')

if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    df = pd.read_table(STATES_INPUT, header=0, index_col=0)
    df.index = df.index.map(str)
    date_column = 'Year'

    country2loc = {'China': 'Asia', 'South Korea': 'Asia',
                   'Myanmar': 'Asia', 'Philippines': 'Asia', 'Thailand': 'Asia',
                   'India': 'Indian subcontinent', 'Nepal': 'Indian subcontinent', 'Pakistan': 'Indian subcontinent',
                   #
                   'Israel': 'Europe', 'Yemen': 'Horn of Africa', 'Georgia': 'Europe',
                   #
                   'Austria': 'Europe', 'Belgium': 'Europe', 'Cyprus': 'Europe', 'Denmark': 'Europe',
                   'Finland': 'Europe', 'France': 'Europe', 'Germany': 'Europe', 'Greece': 'Europe',
                   'Italy': 'Europe', 'Luxembourg': 'Europe', 'Netherlands': 'Europe',
                   'Norway': 'Europe', 'Portugal': 'Europe', 'Spain': 'Europe', 'Sweden': 'Europe',
                   'Switzerland': 'Europe', 'United Kingdom': 'Europe',
                   'Czech Rep.': 'Europe', 'Slovakia': 'Europe', 'Poland': 'Europe', 'Romania': 'Europe',
                   #
                   'Russian Federation': 'Europe', 'Ukraine': 'Europe',
                   #
                   'Australia': 'Australia',
                   #
                   'United States': 'North America',
                   #
                   'Cuba': 'Central America', 'Honduras': 'Central America',
                   #
                   'Argentina': 'South America', 'Brazil': 'South America', 'Uruguay': 'South America',
                   'Venezuela': 'South America',
                   #
                   'Central Africa': 'Central Africa', 'Dem. Rep. of Congo': 'Central Africa',
                   'Zambia': 'Central Africa',
                   #
                   'Burundi': 'East Africa', 'Kenya': 'East Africa', 'Rwanda': 'East Africa', 'Tanzania': 'East Africa',
                   'Uganda': 'East Africa',
                   #
                   'Djibouti': 'Horn of Africa', 'Eritrea': 'Horn of Africa', 'Ethiopia': 'Horn of Africa',
                   'Somalia': 'Horn of Africa', 'Sudan': 'Horn of Africa',
                   #
                   'South Africa': 'South Africa', 'Swaziland': 'South Africa',
                   #
                   'Botswana': 'Southern Africa ex SA', 'Malawi': 'Southern Africa ex SA',
                   'Mozambique': 'Southern Africa ex SA', 'Zimbabwe': 'Southern Africa ex SA',
                   #
                   'Cameroon': 'West Africa', 'Equatorial Guinea': 'West Africa', 'Gabon': 'West Africa',
                   'Mali': 'West Africa', 'Nigeria': 'West Africa', 'Senegal': 'West Africa'}

    iso32loc = {Country.get_iso3_country_code_fuzzy(c)[0]: loc for (c, loc) in country2loc.items()
                if Country.get_iso3_country_code_fuzzy(c)[0]}


    def get_location(_, sa=False):
        if pd.isna(_):
            return None
        loc = iso32loc[_]
        if not sa:
            return loc
        if 'Indian' in loc:
            return 'Asia'
        if 'Europe' in loc or 'Asia' in loc:
            return loc
        return Country.get_country_info_from_iso3(_)['Country or Area']


    df['Loc_SA'] = df['Country ISO3'].map(lambda _: get_location(_, True))
    df['Loc'] = df['Country ISO3'].map(get_location)
    df.to_csv(STATES_INPUT_LOC, sep='\t')

    generate_map(data=STATES_INPUT_LOC, country='Country', location='Loc', tree=TREE_NWK,
                 html=os.path.join(DATA_DIR, 'maps', 'geo_map.html'))

    model = F81
    prediction_method = MARGINAL_APPROXIMATION
    pastml_pipeline(data=STATES_INPUT_LOC, tree=TREE_NWK,
                    html_compressed=os.path.join(DATA_DIR, 'maps', 'map_Loc.html'),
                    model=model, verbose=False, columns=['Loc'],
                    prediction_method=prediction_method, work_dir=os.path.join(DATA_DIR, 'pastml'), date_column='Year')

    mutation = 'RT:M184V'

    pastml_pipeline(data=STATES_INPUT_LOC, tree=TREE_NWK,
                    html_compressed=os.path.join(DATA_DIR, 'maps', mutation, 'map_{}.html'.format(mutation)),
                    model=model, verbose=False, columns=[mutation], prediction_method=prediction_method,
                    work_dir=os.path.join(DATA_DIR, 'pastml'), date_column='Year')

    pastml_pipeline(data=STATES_INPUT_LOC, tree=TREE_NWK,
                    html_compressed=os.path.join(DATA_DIR, 'maps', 'map_Loc_{}.html'.format(mutation)),
                    model=model, verbose=False, columns=[mutation, 'Loc'],
                    prediction_method=prediction_method,
                    name_column='Loc',
                    work_dir=os.path.join(DATA_DIR, 'pastml'), date_column='Year')

    tree = Tree(TREE_NWK, format=3)
    for _ in tree.traverse():
        if _.name == 'node_3596':
            _.write(outfile=SA_SUBTREE_NWK, format=3)

    generate_map(data=STATES_INPUT_LOC, country='Country', location='Loc_SA', tree=SA_SUBTREE_NWK,
                 html=os.path.join(DATA_DIR, 'maps', 'geo_map_SA.html'))

    pastml_pipeline(data=STATES_INPUT_LOC, tree=SA_SUBTREE_NWK,
                    html_compressed=os.path.join(DATA_DIR, 'maps',
                                                 'map_Loc_{}_SA.html'.format(mutation)),
                    model=model, verbose=False, columns=[mutation, 'Loc_SA'],
                    prediction_method=prediction_method,
                    name_column='Loc_SA',
                    work_dir=os.path.join(DATA_DIR, 'pastml'), date_column='Year',
                    column2parameters={'RT:M184V': {'scaling factor': 4.27143674, 'epsilon': 3e-03,
                                                    "sensitive": 0.93371409, "resistant": 0.06628591},
                                       # 'Loc_SA': {'scaling factor': 6.49557354, 'epsilon': 3e-03,
                                       #            "Ethiopia": 0.03704372,
                                       #            "Djibouti": 0.00645968, " Brazil": 0.01283314, " Europe": 0.50565032,
                                       #            " Uganda": 0.05331524, " Senegal": 0.03236439,
                                       #            " Democratic Republic of the Congo": 0.02581537,
                                       #            " Burundi": 0.08141071, " Sudan": 0.02602897, " Eritrea": 0.00651429,
                                       #            " Uruguay": 0.00646014, " Argentina": 0.03883105,
                                       #            " Yemen": 0.01936169, " Kenya": 0.01290804, " Cuba": 0.01980251,
                                       #            " Asia": 0.01290804, " South Africa": 0.00651588,
                                       #            " United Republic of Tanzania": 0.05053965,
                                       #            " United States of America": 0.01936169, " Zambia": 0.02587549}
                                       },
                    tip_size_threshold=10
                    )

    for year in (2011, 2006, 2001, 1996, 1991, 1986):
        tree = remove_certain_leaves(tree, to_remove=lambda node: pd.isnull(df.loc[node.name, date_column])
                                                                  or df.loc[node.name, date_column] > year)
        nwk = YEAR_SUBTREE_NWK.format(year)
        tree.write(outfile=nwk, format=3)

        pastml_pipeline(data=STATES_INPUT_LOC, tree=nwk,
                        html_compressed=os.path.join(DATA_DIR, 'maps', mutation,
                                                     'map_{}_{}.html'.format(mutation, year)),
                        model=model, verbose=False, columns=[mutation],
                        prediction_method=prediction_method,
                        work_dir=os.path.join(DATA_DIR, 'pastml'), date_column='Year')
