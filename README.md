# cytopast

__cytopast__ is a python3 module that creates zoomable HTML visualisation of phylogenetic trees with annotated nodes.

Given a tree and its node annotations, it can either visualise them as-is, 
or apply [PASTML](https://github.com/saishikawa/PASTML) to infer ancestral states based on the tip states. 

The states are visualised as different colours of the tree nodes using [Cytoscape.js](http://js.cytoscape.org/)

# Input data
As an input, one needs to provide a **rooted** phylogenetical tree in [newick](https://en.wikipedia.org/wiki/Newick_format) format,
and a table containing tip states, 
in tab-delimited (by default) or csv format (to be specified with *--data_sep ,* option).

### Example
Let's assume that the tree and annotation files are in the Downloads folder, 
and are named respectively tree.nwk and states.csv.

The states.csv is a comma-separated file, containing tip ids in the first column, 
and several named columns, including *Location*, i.e.:


Tip_id | ... | Location | ...
----- |  ----- | ----- | -----
1 | ... | Africa | ...
2 | ... | Asia | ...
3 | ... | Africa | ...
... | ... | ... | ...


# Try it online
Try it at [pastml.pasteur.fr](https://pastml.pasteur.fr)

# Run it on your computer

There are 2 alternative ways to run cytopast on your computer: with [docker](https://hub.docker.com/), or in python3.

## Run with docker

### Basic usage
```bash
docker run -v <path_to_the_folder_containing_the_tree_and_the_annotations>:/data:rw -t evolbioinfo/pastml --tree /data/<tree_file> --data /data/<annotation_file> --columns <one_or_more_column_names> --html_compressed /data/<map_name>
```

For example, to reconstruct and visualise the ancestral Location states, 
one needs to run the following command:

```bash
docker run -v ~/Downloads:/data:rw -t evolbioinfo/pastml --tree /data/tree.nwk --data /data/states.csv --data_sep , --columns Location --html_compressed /data/location_map.html
```

This will produce a file location_map.html in the Downloads folder, 
that can be viewed with a browser.


### Help

To see advanced options, run
```bash
docker run -t evolbioinfo/pastml -h
```

## Run in python3

We strongly recommend installing cytopast for python via [conda](https://conda.io/docs/), following the procedure described below:

### Installing with conda

Once you have conda installed create an environment for cytopast with python3, gcc and gsl:

```bash
conda create --name cytopast python=3 gcc gsl
```

Then activate it:
```bash
source activate cytopast
```

Then install cytopast in it:

```bash
pip install cytopast
```

### Installing without conda

Install [GNU GSL](https://www.gnu.org/software/gsl/), following the instructions provided on GSL website.

Then install cytopast:

```bash
pip3 install cytopast
```

### Basic usage in a command line
If you installed cytopast via conda, do not forget to first activate the dedicated environment, e.g.

```bash
source activate cytopast
```

To run cytopast:

```bash
cytopast --tree <path/to/tree_file.nwk> --data <path/to/annotation_file.tab> --columns <one_or_more_column_names> --html_compressed <path/to/output/map.html>
```

### Help

To see advanced options, run:
```bash
cytopast -h
```

### Basic usage in python3
```python
from cytopast.pastml_analyser import pastml_pipeline

# Path to the table containing tip/node annotations, in csv or tab format
data = "/path/to/the/table/eg/data.csv"

# Path to the tree in newick format
tree = "/path/to/the/tree/eg/tree.nwk"

# Columns present in the annotation table,
# for which we want to reconstruct ancestral states
columns = ['Location', 'Resistant_or_not']

# Path to the output compressed map visualisation
html_compressed = "/path/to/the/future/map/eg/map.html"

# Path to the output tree visualisation
html = "/path/to/the/future/tree/visualisation/eg/tree.html"

pastml_pipeline(data=data, data_sep=',', columns=columns, name_column='Location',
                tree=tree,
                html_compressed=html_compressed, html=html, 
                verbose=True)
```

### Examples

See the [examples folder](https://github.com/evolbioinfo/cytopast/tree/master/examples) for ideas :)
