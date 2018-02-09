# cytopast

__cytopast__ is a python3 module that creates zoomable HTML visualisation of phylogenetic trees with annotated nodes.

Given a tree and its node annotations, it can either visualise them as-is, 
or apply [PASTML](https://github.com/saishikawa/PASTML) to infer ancestral states based on the tip states. 

The states are visualised as different colours of the tree nodes using [Cytoscape.js](http://js.cytoscape.org/)

## Installation
```bash
pip3 install cytopast
```

## Basic usage in python3
```python
from cytopast.pastml_analyser import pastml_pipeline

# Path to the table containing tip/node annotations, in csv or tab format
data = "/path/to/the/table/eg/data.csv"

# Path to the tree in newick format
tree = "/path/to/the/tree/eg/tree.nwk"

# Columns present in the annotation table,
# for which we want to reconstruct ancestral states
columns = ['Location', 'Resistant_or_not']

# Columns present in the annotation table,
# for which we want to copy existing annotations from the annotation table,
# without inferring ancestral states
copy_columns = ['Sex']

# Path to the output compressed map visualisation
html_compressed = "/path/to/the/future/map/eg/map.html"

# Path to the output tree visualisation
html = "/path/to/the/future/tree/visualisation/eg/tree.html"

pastml_pipeline(data=data, data_sep=',', columns=columns, name_column='Location',
                tree=tree,
                html_compressed=html_compressed, html=html, 
                verbose=True)
```

## Basic usage from console
```bash
cytopast --tree /path/to/the/tree.nwk --data /path/to/the/annotation/data.txt --data_sep , \
--html /path/to/the/output/visualisation/of/the/tree.html \
--html_compressed /path/to/the/output/visualisation/of/the/compressed/map.html \
--columns Location Resistant_or_not --name_column Location --verbose
```

## Options

```
usage: cytopast [-h] -d DATA [-s DATA_SEP] [-i ID_INDEX]
                [-c [COLUMNS [COLUMNS ...]]]
                [--copy_columns [COPY_COLUMNS [COPY_COLUMNS ...]]] -t TREE
                [-m MODEL] [--work_dir WORK_DIR] [--cache] [-n NAME_COLUMN]
                [-a] [-o OUT_DATA] [-p HTML_COMPRESSED] [-l HTML] [-v]

Visualisation of annotated phylogenetic trees (as html maps).

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         print information on the progress of the analysis

annotation-related arguments:
  -d DATA, --data DATA  the annotation file in tab/csv format with the first
                        row containing the column names.
  -s DATA_SEP, --data_sep DATA_SEP
                        the column separator for the data table. By default is
                        set to tab, i.e. for tab file. Set it to ',' if your
                        file is csv.
  -i ID_INDEX, --id_index ID_INDEX
                        the index of the column in the data table that
                        contains the tree tip names, indices start from zero
                        (by default is set to 0).
  -c [COLUMNS [COLUMNS ...]], --columns [COLUMNS [COLUMNS ...]]
                        names of the data table columns that contain states to
                        be analysed with PASTML. If neither columns nor
                        copy_columns are specified, then all columns will be
                        considered for PASTMl analysis.
  --copy_columns [COPY_COLUMNS [COPY_COLUMNS ...]]
                        names of the data table columns that contain states to
                        be copied as-is, without applying PASTML (the missing
                        states will stay unresolved).

tree-related arguments:
  -t TREE, --tree TREE  the input tree in newick format.

ancestral-state inference-related arguments:
  -m MODEL, --model MODEL
                        the evolutionary model to be used by PASTML (can be JC
                        or F81).
  --work_dir WORK_DIR   the working dir for PASTML to put intermediate files
                        into (if not specified a temporary dir will be
                        created).
  --cache               if set, the results of previous PASTML runs on this
                        data will be reused when possible

visualisation-related arguments:
  -n NAME_COLUMN, --name_column NAME_COLUMN
                        name of the data table column to be used for node
                        names in the compressed map visualisation(must be one
                        of those specified in columns or copy_columns if they
                        are specified).If the data table contains only one
                        column it will be used by default.
  -a, --all             Keep all the nodes in the compressed map
                        visualisation, even the minor ones.

output-related arguments:
  -o OUT_DATA, --out_data OUT_DATA
                        the output annotation file with the states inferred by
                        PASTML.
  -p HTML_COMPRESSED, --html_compressed HTML_COMPRESSED
                        the output summary map visualisation file (html).
  -l HTML, --html HTML  the output tree visualisation file (html).
```