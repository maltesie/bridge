# Bridge
PyMol plugin for hydrogen bond and water wire analysis of crystal structures and MD simulations.

## Copyright Notice and Disclaimer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Author: Malte Siemers, Freie Universität Berlin 
   
If you use this software or anything it produces for work to be published,
please cite:
   
Malte Siemers, Michalis Lazaratos, Konstantina Karathanou, Federico Guerra, Leonid Brown, and Ana-Nicoleta Bondar.
Bridge: A graph-based algorithm to analyze dynamic H-bond networks in membrane proteins, Journal of Chemical Theory and Computation, 2019.

and

Federico Guerra, Malte Siemers, Christopher Mielack, and Ana-Nicoleta Bondar<br/>
Dynamics of Long-Distance Hydrogen-Bond Networks in Photosystem II<br/>
The Journal of Physical Chemistry B 2018 122 (17), 4625-4641 <br/>

## Installation Guide

Bridge comes as a pymol plugin and depends on MDAnalysis (0.19.2+). 
To install pymol, we recommend using the miniconda python distribution 
(https://docs.conda.io/en/latest/miniconda.html). If miniconda is set up 
correctly you can install pymol and MDAnalysis with:

> conda install -c schrodinger pymol <br/>
> conda install -c conda-forge mdanalysis

To open pymols plugin manager, click

> Plugin->Plugin Manager

In the plugin manager, go to the 'Install New Plugin' tab and click 
'Choose file ...'. Select the .zip file containing Bridge and restart
pymol. Now start Bridge by clicking

> Plugin->Bridge
