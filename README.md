# Bridge
PyMol plugin for hydrogen bond and water wire analysis of crystal structures and MD simulations.

## State of Publication
The code for this project will be publically available when the article

Malte Siemers, Michalis Lazaratos, Konstantina Karathanou,
Federico Guerra, Leonid Brown, and Ana-Nicoleta Bondar. 
Bridge: A graph-based algorithm to analyze dynamic H-bond networks 
in membrane proteins, Journal of Chemical Theory and Computation, reviewed August 2019.

is published.

## Installation Guide

Bridge comes as a pymol plugin and depends on MDAnalysis (0.19.2+). 
To install pymol, we recommend using the miniconda python distribution 
(https://docs.conda.io/en/latest/miniconda.html). If miniconda is set up 
correctly you can install pymol and  MDAnalysis with:

conda install -c schrodinger pymol
conda install -c conda-forge mdanalysis

To open pymols plugin manager, click

Plugin->Plugin Manager

In the plugin manager, go to the 'Install New Plugin' tab and click 
'Choose file ...'. Select the .zip file containing Bridge and restart
pymol. Now start Bridge by clicking

Plugin->Bridge
