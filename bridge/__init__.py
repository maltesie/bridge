#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#    Author: Malte Siemers, Freie Universität Berlin 
#   
#    If you use this software or anything it produces for work to be published,
#    please cite:
#    
#    Malte Siemers, Michalis Lazaratos, Konstantina Karathanou,
#    Federico Guerra, Leonid Brown, and Ana-Nicoleta Bondar. 
#    Bridge: A graph-based algorithm to analyze dynamic H-bond networks in
#    membrane proteins, Journal of Chemical Theory and Computation 15 (12) 6781-6798
#
#    and
#
#    Federico Guerra, Malte Siemers, Christopher Mielack, and Ana-Nicoleta Bondar
#    Dynamics of Long-Distance Hydrogen-Bond Networks in Photosystem II
#    The Journal of Physical Chemistry B 2018 122 (17), 4625-4641


from __future__ import absolute_import
from __future__ import print_function

import os
wa, hba = None, None
add_hba, add_wa = [], []
static_canvas, current_figure, current_data = None, None, ''
backbone_backbone = False

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Bridge', run_plugin_gui)

def run_plugin_gui():
    '''
    Open our custom dialog
    '''
    # entry point to PyMOL's API
    #from pymol import cmd
    from .mdhbond import HbondAnalysis
    from .mdhbond import WireAnalysis
    from .mdhbond.helpfunctions import donor_names_global, acceptor_names_global

    # pymol.Qt provides the PyQt5 interface, but may support PyQt4
    # and/or PySide as well
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi
    from pymol.Qt.utils import getSaveFileNameWithExt
    import numpy as np
    from matplotlib.figure import Figure
    import matplotlib.pyplot as plt

    from matplotlib.backends.qt_compat import is_pyqt5
    if is_pyqt5():
        from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
    else:
        from matplotlib.backends.backend_qt4agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar) 
    
    class ApplicationWindow(QtWidgets.QMainWindow):
        def __init__(self):
            super().__init__()
            self._main = QtWidgets.QWidget()
            self.setCentralWidget(self._main)
            self.layout = QtWidgets.QVBoxLayout(self._main)
    
            self.static_canvas = FigureCanvas(Figure(figsize=(20, 15)))
            self.layout.addWidget(self.static_canvas)
            self.addToolBar(NavigationToolbar(self.static_canvas, self))
    
    def getOpenFileNameWithExt(*args, **kwargs):
        """
        Return a file name, append extension from filter if no extension provided.
        """
        import os, re
    
        fname, filter = QtWidgets.QFileDialog.getOpenFileName(*args, **kwargs)
    
        if not fname:
            return ''
    
        if '.' not in os.path.split(fname)[-1]:
            m = re.search(r'\*(\.[\w\.]+)', filter)
            if m:
                # append first extension from filter
                fname += m.group(1)
    
        return fname
    
    def getOpenFileNamesWithExt(*args, **kwargs):
        """
        Return a file name, append extension from filter if no extension provided.
        """
        import os, re
    
        fnames, filter = QtWidgets.QFileDialog.getOpenFileNames(*args, **kwargs)
    
        if not fnames:
            return ''
        
        fnames_out = []
        for fname in fnames:
            if '.' not in os.path.split(fname)[-1]:
                m = re.search(r'\*(\.[\w\.]+)', filter)
                if m:
                    # append first extension from filter
                    fname += m.group(1)
            fnames_out.append(fname)
        
        return fnames_out
    
    def make_legend_list(line):
        if line == '': return None
        legend = line.split(',')
        for i in range(len(legend)):
            legend[i] = legend[i].strip(' ')
        return legend
    
    # create a new Window
    dialog = QtWidgets.QMainWindow()
    
    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'mainwidget.ui')
    form = loadUi(uifile, dialog)
    dialog.setCentralWidget(form.scrollArea)
    
    #graph_app = ApplicationWindow()
    
    uifile_plotting = os.path.join(os.path.dirname(__file__), 'plotting_unified.ui')
    dialog_plotting = QtWidgets.QDialog()
    form_plotting = loadUi(uifile_plotting, dialog_plotting)
    
    uifile_graph_bonds = os.path.join(os.path.dirname(__file__), 'advanced_graph_visualization_bonds.ui')
    dialog_graph_bonds = QtWidgets.QDialog()
    form_graph_bonds = loadUi(uifile_graph_bonds, dialog_graph_bonds)
    
    uifile_graph_wires = os.path.join(os.path.dirname(__file__), 'advanced_graph_visualization_wires.ui')
    dialog_graph_wires = QtWidgets.QDialog()
    form_graph_wires = loadUi(uifile_graph_wires, dialog_graph_wires)
    
    global static_canvas
    static_canvas = FigureCanvas(Figure(figsize=(6.5, 4)))
    form_plotting.gridLayout_22.addWidget(static_canvas)
    
    uifile_error = os.path.join(os.path.dirname(__file__), 'error_line_ok.ui')
    dialog_error = QtWidgets.QDialog()
    form_error = loadUi(uifile_error, dialog_error)
    
    def write_license():
        print("\n\n    Bridge\n\n    This program is free software: you can redistribute it and/or modify\n    it under the terms of the GNU General Public License as published by\n    the Free Software Foundation, either version 3 of the License, or\n    (at your option) any later version.\n\n    This program is distributed in the hope that it will be useful,\n    but WITHOUT ANY WARRANTY; without even the implied warranty of\n    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n    GNU General Public License for more details.\n\n    You should have received a copy of the GNU General Public License\n    along with this program.  If not, see <https://www.gnu.org/licenses/>.\n\n    Author: Malte Siemers, Freie Universität Berlin\n\n    If you use this software or anything it produces for work to be published,\n    please cite:\n\n    Malte Siemers, Michalis Lazaratos, Konstantina Karathanou,\n    Federico Guerra, Leonid Brown, and Ana-Nicoleta Bondar.\n    Bridge: A graph-based algorithm to analyze dynamic H-bond networks\n    in membrane proteins,\n    Journal of Chemical Theory and Computation 2019, 15 (12) 6781-6798\n\n    and\n\n    Federico Guerra, Malte Siemers, Christopher Mielack, and Ana-Nicoleta Bondar\n    Dynamics of Long-Distance Hydrogen-Bond Networks in Photosystem II\n    The Journal of Physical Chemistry B 2018 122 (17), 4625-4641")
    
    def draw_figure_to_canvas(f, print_string):  
        global current_figure
        global current_data
        current_data = print_string
        w,h = f.get_size_inches()
        try: scale = float(form_plotting.line_scale.text())
        except: scale = 1.0
        f.set_size_inches(scale*w,scale*h)
        current_figure = FigureCanvas(f)
        f.set_canvas(static_canvas)
        f.draw(f.canvas.get_renderer())
        static_canvas.update()
        if not dialog.isVisible():
            dialog.show()
    
    def save_plot():
        filename = getSaveFileNameWithExt(dialog_plotting, 'Save', filter='Portable network graphic (*.png);;Vector graphic (*.eps)')
        if filename == '': return
        end = filename.split('.')[-1]
        if end == 'eps': plt.text.usetex = True
        current_figure.print_figure(filename, format=end, dpi=300)
    
    def save_data():
        filename = getSaveFileNameWithExt(dialog_plotting, 'Save', filter='ASCII data file (*.txt)')
        with open(filename, 'w') as file:
            file.write(current_data)
    
    def error(text):
        form_error.label_error.setText(text)
        load_error()
    
    def load_error():
        dialog_error.show()
    
    def close_error():
        dialog_error.close()
        
    form_error.button_ok.clicked.connect(close_error)
    
    def enable_tabs():
        if hba is not None:
            form_plotting.tab.setEnabled(True)
        if wa is not None:
            form_plotting.tab_2.setEnabled(True)
    
    def save_wires_data():
        filename = getSaveFileNameWithExt(dialog_plotting, 'Save', filter='ASCII data file (*.txt)')
        if not filename: return
        save_str = 'wire_partner_1 wire_partner_2 occupancy\n\n'
        for bond in wa.filtered_results:
            ba, bb = bond.split(':')
            save_str += ba + ' ' + bb + ' ' + str(np.round(wa.filtered_results[bond].mean(), 3))+'\n'
        with open(filename, 'w') as file:
            file.write(save_str)
            
    def save_wires_data_advanced():
        filename = getSaveFileNameWithExt(dialog_plotting, 'Save', filter='ASCII data file (*.txt)')
        if not filename: return
        try:
            weight = float(form_graph_wires.lineEdit_wires_color_weight.text())
        except: 
            weight = None
        if form_graph_wires.radioButton_bonds_colors_segments.isChecked():
            save_str = 'wire_partner_1 wire_partner_2 occupancy average_water\n\n'
            for bond in wa.filtered_results:
                ba, bb = bond.split(':')
                save_str += ba + ' ' + bb + ' ' + str(np.round(wa.filtered_results[bond].mean(), 2)) + ' ' + str(np.round(wa.wire_lengths[bond][wa.wire_lengths[bond]!=np.inf].mean(), 2)) +'\n'
        elif form_graph_wires.radioButton_bonds_betweenness.isChecked():
            save_str = 'node betweenness_centrality\n\n'
            betweenness = wa.compute_centrality(centrality_type='betweenness', weight=weight)
            for node in betweenness:
                save_str += node + ' ' + str(np.round(betweenness[node], 2))+'\n'
        elif form_graph_wires.radioButton_bonds_degree.isChecked():
            save_str = 'node degree_centrality\n\n'
            degree = wa.compute_centrality(centrality_type='degree', weight=weight)
            for node in degree:
                save_str += node + ' ' + str(np.round(degree[node], 2))+'\n'
        with open(filename, 'w') as file:
            file.write(save_str)
    
    def save_bonds_data():
        filename = getSaveFileNameWithExt(dialog_plotting, 'Save', filter='ASCII data file (*.txt)')
        if not filename: return
        save_str = 'H_bond_partner_1 H_bond_partner_2 occupancy\n\n'
        for bond in hba.filtered_results:
            ba, bb = bond.split(':')
            save_str += ba + ' ' + bb + ' ' + str(np.round(hba.filtered_results[bond].mean(), 3))+'\n'
        with open(filename, 'w') as file:
            file.write(save_str)
            
    def save_bonds_data_advanced():
        filename = getSaveFileNameWithExt(dialog_plotting, 'Save', filter='ASCII data file (*.txt)')
        if not filename: return
        try:
            weight = float(form_graph_bonds.lineEdit_bonds_color_weight.text())
        except: 
            weight = None
        if form_graph_bonds.radioButton_bonds_colors_segments.isChecked():
            save_str = 'H_bond_partner_1 H_bond_partner_2 occupancy\n\n'
            for bond in hba.filtered_results:
                ba, bb = bond.split(':')
                save_str += ba + ' ' + bb + ' ' + str(np.round(hba.filtered_results[bond].mean(), 2)) +'\n'
        elif form_graph_bonds.radioButton_bonds_betweenness.isChecked():
            save_str = 'node betweenness_centrality\n\n'
            betweenness = hba.compute_centrality(centrality_type='betweenness', weight=weight)
            for node in betweenness:
                save_str += node + ' ' + str(np.round(betweenness[node], 2))+'\n'
        elif form_graph_bonds.radioButton_bonds_degree.isChecked():
            save_str = 'node degree_centrality\n\n'
            degree = hba.compute_centrality(centrality_type='degree', weight=weight)
            for node in degree:
                save_str += node + ' ' + str(np.round(degree[node], 2))+'\n'
        with open(filename, 'w') as file:
            file.write(save_str)
    
    def load_advanced_wires():
        dialog_graph_wires.show()
        
    def load_advanced_bonds():
        dialog_graph_bonds.show() 
    
    
    def load_plotting_wires():
        enable_tabs()
        form_plotting.tabWidget.setCurrentWidget(form_plotting.tab_2)
        dialog_plotting.show()
        
    def load_plotting_bonds():
        enable_tabs()
        form_plotting.tabWidget.setCurrentWidget(form_plotting.tab)
        dialog_plotting.show()
    
    def toggle_bonds_angle():
        form.line_bonds_angle.setEnabled(form.checkBox_angle.isChecked())
    
    def toggle_wire_angle():
        form.line_wire_angle.setEnabled(form.checkBox_wires_angle.isChecked())
    
    def get_applied_filter_hba():
        exp = hba.applied_filters['single_path']
        occ = hba.applied_filters['occupancy']
        paths = hba.applied_filters['shortest_paths']
        comp = hba.applied_filters['connected_component']
        segm = hba.applied_filters['segnames']
        resn = hba.applied_filters['resnames']
        hyd = hba.applied_filters['shells']
        back = hba.applied_filters['backbone']
        return exp, occ, comp, paths, segm, resn, hyd, back
    
    def get_applied_filter_wa():
        exp = wa.applied_filters['single_path']
        occ = wa.applied_filters['occupancy']
        paths = wa.applied_filters['shortest_paths']
        comp = wa.applied_filters['connected_component']
        segm = wa.applied_filters['segnames']
        resn = wa.applied_filters['resnames']
        avg_least = wa.applied_filters['avg_least_bonds']
        return exp, occ, comp, paths, segm, resn, avg_least
    
    def browse_filename_structure():
        filename = getOpenFileNameWithExt(
            dialog, 'Open', filter='PSF File (*.psf);;PDB File (*.pdb);;All Files (*.*)')
        if filename:
            form.line_bonds_structure.setText(filename)
    
    def browse_filename_trajectories():
        filename = getOpenFileNamesWithExt(
            dialog, 'Open', filter='DCD File (*.dcd);;All Files (*.*)')
        if filename:
            fnames = ', '.join(filename)
            form.line_bonds_trajectories.setText(fnames)
    
    
    def init_bonds_analysis():
        global hba
        form.button_bonds_init.setText('Working...')
        form.button_bonds_init.repaint()
        ad = [donor for donor in form.line_bonds_donors.text().replace(' ','').split(',') if donor != '']
        aa = [acceptor for acceptor in form.line_bonds_acceptors.text().replace(' ','').split(',') if acceptor != '']
        ions = [ion for ion in form.line_bonds_ions.text().replace(' ','').split(',') if ion != '']
        form.group_bonds_filter.setEnabled(False)
        form.group_bonds_filter.repaint()
        #form.group_bonds_compute.setEnabled(True)
        #form.group_bonds_draw.setEnabled(True)
        form.groupBox_2D_bonds.setEnabled(False)
        form.groupBox_2D_bonds.repaint()
        form.groupBox_plotting_bonds.setEnabled(False)
        form.groupBox_plotting_bonds.repaint()
        #form.group_bonds_save.setEnabled(True)
        form.button_bonds_save.setEnabled(False)
        form.button_bonds_save.repaint()
        form.group_set_bonds.setEnabled(False)
        form.group_set_bonds.repaint()
        form.textedit_bonds_active_filters.clear()
        form.repaint()
        set_params_to_default_bonds()
        try:
            step = form.line_bonds_step.text()
            stop = form.line_bonds_stop.text()
            start = form.line_bonds_start.text()
            if step == '': step=1
            else: step = int(step)
            if stop == '': stop=None
            else: stop = int(stop)
            if start == '': start=None
            else: start = int(start)
            trajectory = form.line_bonds_trajectories.text().split(', ')
            if trajectory == ['']: trajectory = None
            donors_without = form.checkBox_bonds_donors_without_hydrogen.isChecked()
            check_angle=True
            if donors_without: 
                form.checkBox_angle.setChecked(False)
                check_angle=False
            try:
                hba = HbondAnalysis(selection=form.line_bonds_selection.text(), structure=form.line_bonds_structure.text(),
                                trajectories=trajectory, check_angle=check_angle, residuewise=form.checkbox_bonds_residuewise.isChecked(), 
                                additional_donors=ad, additional_acceptors=aa, ions=ions, add_donors_without_hydrogen=donors_without)
                try: hba.add_missing_residues = int(form.line_bonds_missing_residues.text())
                except: error('Missing residues number not set. Cant convert input to integer.')
            except AssertionError:
                error('Couldnt initialize analysis. Are there hydrogens in the structure?')
                form.button_bonds_init.setText('Initialize')
                form.button_bonds_init.repaint()
                return
            form.group_set_bonds.setEnabled(True)
            if not form.checkbox_bonds_residuewise.isChecked(): 
                form_plotting.groupBox.setEnabled(True)
                form.groupBox_filter_backbone.setEnabled(True)
                form_plotting.checkbox_bonds_heatmap_residuewise.setEnabled(True)
            else: 
                form_plotting.groupBox.setEnabled(False)
                form.groupBox_filter_backbone.setEnabled(False)
                form_plotting.checkbox_bonds_heatmap_residuewise.setEnabled(False)
                form_plotting.checkbox_bonds_heatmap_residuewise.setChecked(False)
                toggle_hbonds_heatmap_lines()
            form_plotting.groupBox.repaint()
            form.button_bonds_init.setText('Initialize')
            form.button_bonds_init.repaint()
        except:
            form.button_bonds_init.setText('Initialize')
            form.button_bonds_init.repaint()
            error("One or more input fields in the Initialize area are of unexpected type or empty.")
        
    def set_bonds():
        global backbone_backbone
        try:
            step = form.line_bonds_step.text()
            stop = form.line_bonds_stop.text()
            start = form.line_bonds_start.text()
            if step == '': step=1
            else: step = int(step)
            if stop == '': stop=None
            else: stop = int(stop)
            if start == '': start=None
            else: start = int(start)
            dist = float(form.line_bonds_dist.text())
            angle = float(form.line_bonds_angle.text())
            check_angle = form.checkBox_angle.isChecked()
            if form.checkBox_bonds_donors_without_hydrogen.isChecked() and check_angle: 
                error('Donors without hydrogen can only be used with no angle check! Proceeding with angle check set to false')
                form.checkBox_angle.setChecked(False)
                check_angle = False
            ad = [donor if donor!='' else list(donor_names_global)[0] for donor in form.line_bonds_donors.text().replace(' ','').split(',')]
            aa = [acceptor if acceptor!='' else list(acceptor_names_global)[0] for acceptor in form.line_bonds_acceptors.text().replace(' ','').split(',')]
            hba.check_angle = check_angle
            hba.distance = dist
            hba.cut_angle = angle
            hba._trajectory_slice = slice(start if isinstance(start, int) else None, stop if isinstance(stop, int) else None, step)
            hba.nb_frames = len([0 for i in hba._universe.trajectory[hba._trajectory_slice]])
            hba.donor_names = hba.donor_names.union(ad)
            hba.acceptor_names = hba.acceptor_names.union(aa)
            backbone_backbone = form.checkBox_bonds_exclude_backbone.isChecked()
            hba._exclude_backbone_backbone = backbone_backbone
        except:
            error('One or more input fields in the Search area are of unexpected type.')
        form.bonds_set.setText('Working...')
        form.bonds_set.repaint()
        if form.radio_in_selection.isChecked():
            hba.set_hbonds_in_selection(exclude_backbone_backbone=backbone_backbone)
        elif form.radio_around.isChecked():
            try:
                try: around_rad = float(form.line_around_value.text())
                except: error('Wrong input type for the around value')
                hba.set_hbonds_in_selection_and_water_around(around_rad, exclude_backbone_backbone=backbone_backbone)
            except:
                error("Couldt run the algorithm. Is water present?")
                form.bonds_set.setText('Set')
                form.bonds_set.repaint()
                return
        elif form.radio_in_hull.isChecked():
            try: hba.set_hbonds_in_selection_and_water_in_convex_hull(exclude_backbone_backbone=backbone_backbone)
            except: error('Could not run algorithm. Is water present?')
        else:
            hba.set_hbonds_only_water_in_convex_hull()
        clear_all_filters_bonds()
        form.bonds_set.setText('Search')
        form.bonds_set.repaint()
        form.group_bonds_filter.setEnabled(True)
        #form.group_bonds_compute.setEnabled(True)
        #form.group_bonds_draw.setEnabled(True)
        form.groupBox_2D_bonds.setEnabled(True)
        form.groupBox_plotting_bonds.setEnabled(True)
        #form.group_bonds_save.setEnabled(True)
        form.button_bonds_save.setEnabled(True)
    
    def clear_all_filters_bonds():
        form.textedit_bonds_active_filters.clear()
        form.line_bonds_occupancy.setText('')
        form.textedit_bonds_path.setText('')
        form.line_bonds_connected_root.setText('')
        form.line_bonds_path_root.setText('')
        form.line_bonds_path_goal.setText('')
        form.line_bonds_filter_sega.setText('')
        form.line_bonds_filter_segb.setText('')
        form.line_bonds_filter_resa.setText('')
        form.line_bonds_filter_resb.setText('')
        form.line_bonds_order_test.setText('')
        form.checkBox_bonds_backbone.setChecked(False)
        filter_bonds()
        form.repaint()
        
    def set_params_to_default_bonds():
        form.line_bonds_start.clear()
        form.line_bonds_stop.clear()
        form.line_bonds_step.clear()
        form.line_bonds_dist.setText('3.5')
        form.line_bonds_angle.setText('60')
        form.repaint()
    
    def filter_bonds():
        occupancy_cut = form.line_bonds_occupancy.text()
        connected_root = form.line_bonds_connected_root.text()
        path_root = form.line_bonds_path_root.text()
        path_goal = form.line_bonds_path_goal.text()
        path_explicit = form.textedit_bonds_path.toPlainText()
        sega = form.line_bonds_filter_sega.text()
        segb = form.line_bonds_filter_segb.text()
        resa = form.line_bonds_filter_resa.text()
        resb = form.line_bonds_filter_resb.text()
        backbone = form.checkBox_bonds_backbone.isChecked()
        form.textedit_bonds_active_filters.clear()
        use_filtered = False
        order = [int(a) for a in form.line_bonds_order_test.text().replace(' ','').split(',') if a!='']
        if len(order) == 0: order=[6,3,1,2,4,5] 
        try:
            order = [int(a) for a in form.line_bonds_order_test.text().replace(' ','').split(',') if a!='']
            if len(order) == 0: order=[6,3,1,2,4,5] 
        except:
            error('order must be integers separated by commas or empty')
            return
        for a in order:
            if path_explicit != '' and a==3:
                try:
                    args = [node for node in path_explicit.replace(' ','').strip('\n').split(',')]
                    hba.filter_single_path(args, use_filtered)
                    form.textedit_bonds_active_filters.append('Explicit path: '+', '.join(args))
                    use_filtered = True
                except: error('Either one or more of the nodes specified are not in the graph, or the nodes dont build a path in that order!')
            if occupancy_cut != '' and a==1:
                try: 
                    occupancy_cut = float(occupancy_cut)
                    hba.filter_occupancy(occupancy_cut, use_filtered)
                    form.textedit_bonds_active_filters.append('Occupancy: '+str(occupancy_cut))
                    use_filtered = True
                except: error('min. value has to be a float!')
            if connected_root != '' and a==2:
                try:
                    hba.filter_connected_component(connected_root, use_filtered)
                    form.textedit_bonds_active_filters.append('Connected component: '+connected_root)
                    use_filtered = True
                except: error('The specified root is not in the graph!')
            if (path_root != '' and path_goal != '') and a==4:
                try:
                    hba.filter_all_paths(path_root, path_goal, use_filtered)
                    form.textedit_bonds_active_filters.append('Paths from '+path_root+' to '+path_goal)
                    use_filtered = True
                except: error('root and/or target are not in the graph!')
            if (sega != '' or resa != '') and a==5:
                try:
                    if sega != '':
                        if segb == '': segb=None
                        hba.filter_between_segnames(sega, segb, use_filtered)
                        if segb is not None: form.textedit_bonds_active_filters.append('Between segments '+sega+' and '+segb)
                        else: form.textedit_bonds_active_filters.append('Between segment '+sega+' and all')
                        use_filtered = True
                    if resa != '':
                        if resb == '': resb=None
                        hba.filter_between_resnames(resa, resb, use_filtered)
                        if resb is not None: form.textedit_bonds_active_filters.append('Between resnames '+resa+' and '+resb)
                        else: form.textedit_bonds_active_filters.append('Between resname '+resa+' and all')
                        use_filtered = True
                except:
                    pass
            if backbone and a==6:
                hba.filter_only_backbone_bonds(use_filtered)
                form.textedit_bonds_active_filters.append('Only backbone-backbone H Bonds')
                use_filtered = True
            
        if len(hba.filtered_results) == 0:
            error('This combination of filters leads to an empty set of results!')
            use_filtered = False
        if not use_filtered: 
            hba.filtered_graph = hba.initial_graph
            hba.filtered_results = hba.initial_results
            form.textedit_bonds_active_filters.clear()
            
    def compute_bonds_joint():
        joint_occupancy = hba.compute_joint_occupancy()
        form.textedit_bonds_results.clear()
        if joint_occupancy.mean()>0.0:
            form.textedit_bonds_results.append('Joint occupancy: '+str(np.round(joint_occupancy.mean(), 2))+', first frame: '+str(np.argmax(joint_occupancy)))
        else:
            form.textedit_bonds_results.append('Joint occupancy: 0.0')
    
    def compute_i4():
        if form.checkbox_bonds_residuewise.isChecked():
            error('This computation is only possible, if residuewise is set false before initialization.')
        else:
            i4_printable = hba.compute_i4_bonds(return_print_table=True)
            form.textedit_bonds_results.clear()
            form.textedit_bonds_results.setText(i4_printable)
        
    def compute_bonds_residence():
        import numpy as _np
        residence_t = hba.compute_residence_times()
        residence_printable = ''
        for key in residence_t:
            residence_printable += key+ '.'*(45-len(key)) + ' | '+ str(_np.round(residence_t[key],2))+'\n'
        form.textedit_bonds_results.clear()
        form.textedit_bonds_results.setText(residence_printable)
    
    def toggle_bonds_occupancy_segments():
        form_plotting.line_bonds_occupancy_segnames.setEnabled(form_plotting.checkbox_bonds_occupancy_segnames.isChecked())
        form_plotting.line_bonds_occupancy_histogram_legend.setEnabled(not form_plotting.checkbox_bonds_occupancy_segnames.isChecked())
        form_plotting.label_9.setEnabled(not form_plotting.checkbox_bonds_occupancy_segnames.isChecked())
        form_plotting.button_bonds_compare_occupancy.setEnabled(not form_plotting.checkbox_bonds_occupancy_segnames.isChecked())
    
    def toggle_bonds_timeseries_segments():
        form_plotting.line_bonds_timeseries_segnames.setEnabled(form_plotting.checkbox_bonds_timeseries_segnames.isChecked())
        form_plotting.line_bonds_sum_legend.setEnabled(not form_plotting.checkbox_bonds_timeseries_segnames.isChecked())
        form_plotting.label_8.setEnabled(not form_plotting.checkbox_bonds_timeseries_segnames.isChecked())
        form_plotting.button_bonds_compare_bonds_timeseries.setEnabled(not form_plotting.checkbox_bonds_timeseries_segnames.isChecked())
    
    def draw_bonds_graph():
        form.button_bonds_draw_graph.setText('Working...')
        form.button_bonds_draw_graph.repaint()
        hba.draw_graph(draw_edge_occupancy=False, draw_labels=True, node_factor=0.6)
        form.button_bonds_draw_graph.setText('Draw')
        form.button_bonds_draw_graph.repaint()
        
    def draw_bonds_graph_advanced():
        color_d={}
        if form_graph_bonds.line_bonds_colors.text() != '': color_d = {key.split(':')[0]:key.split(':')[1] for key in form_graph_bonds.line_bonds_colors.text().replace(' ','').split(',')}
        centrality = None
        weight = None
        if form_graph_bonds.lineEdit_bonds_color_weight.text() != '': 
            try: weight = float(form_graph_bonds.lineEdit_bonds_color_weight.text())
            except: pass
        if form_graph_bonds.radioButton_bonds_betweenness.isChecked():
            centrality = hba.compute_centrality(centrality_type='betweenness', weight=None)
        elif form_graph_bonds.radioButton_bonds_degree.isChecked():
            centrality = hba.compute_centrality(centrality_type='degree', normalize=True, weight=None)
        if form_graph_bonds.checkBox_bonds_weight_centrality.isChecked() and (form_graph_bonds.radioButton_bonds_betweenness.isChecked() or form_graph_bonds.radioButton_bonds_degree.isChecked()):
            ma = np.max(list(centrality.values()))
            for node in centrality:
                centrality[node] /= ma
        labels = form_graph_bonds.checkBox_bonds_graph_labels.isChecked()
        interregion = form_graph_bonds.checkBox_bonds_interregion.isChecked()
        occupancy = form_graph_bonds.checkBox_bonds_occupancy.isChecked()
        form_graph_bonds.pushButton_draw_graph.setText('Working...')
        form_graph_bonds.pushButton_draw_graph.repaint()
        hba.draw_graph(draw_edge_occupancy=occupancy, highlight_interregion=interregion, centrality=centrality, max_centrality=weight, color_dict=color_d, draw_labels=labels, node_factor=0.6)
        form_graph_bonds.close()
        form_graph_bonds.pushButton_draw_graph.setText('Draw')
        form_graph_bonds.pushButton_draw_graph.repaint()
        
    def draw_bonds_graph_compare():
        filename = getOpenFileNameWithExt(dialog_graph_bonds, 'Open', filter='hydrogen bond analysis file (*.hba);;water wire analysis file (*.wwa)')
        if not filename: return
        hba2 = HbondAnalysis(restore_filename=filename)
        mutations=None
        if form_graph_bonds.lineEdit_mutations.text() != '':
            mutations = form_graph_bonds.lineEdit_mutations.text().replace(' ', '').split(',')
        color_d={}
        if form_graph_bonds.line_bonds_colors.text() != '': color_d = {key.split(':')[0]:key.split(':')[1] for key in form_graph_bonds.line_bonds_colors.text().replace(' ','').split(',')}
        centrality = None
        centrality2 = None
        weight = None
        if form_graph_bonds.lineEdit_bonds_color_weight.text() != '': 
            try: weight = float(form_graph_bonds.lineEdit_bonds_color_weight.text())
            except: pass
        if form_graph_bonds.radioButton_bonds_betweenness.isChecked():
            if form_graph_bonds.checkBox_bonds_weight_centrality.isChecked(): centrality = hba.compute_centrality(centrality_type='betweenness', weight=None)
            centrality2 = hba2.compute_centrality(centrality_type='betweenness', weight=None)
        elif form_graph_bonds.radioButton_bonds_degree.isChecked():
            if form_graph_bonds.checkBox_bonds_weight_centrality.isChecked(): centrality = hba.compute_centrality(centrality_type='degree', normalize=True, weight=None)
            centrality2 = hba2.compute_centrality(centrality_type='degree', weight=None)
        if form_graph_bonds.checkBox_bonds_weight_centrality.isChecked() and (form_graph_bonds.radioButton_bonds_betweenness.isChecked() or form_graph_bonds.radioButton_bonds_degree.isChecked()):
            ma = np.max([np.max(list(centrality2.values())), np.max(list(centrality.values()))])
            for node in centrality2:
                centrality2[node] /= ma
        labels = form_graph_bonds.checkBox_bonds_graph_labels.isChecked()
        interregion = form_graph_bonds.checkBox_bonds_interregion.isChecked()
        occupancy = form_graph_bonds.checkBox_bonds_occupancy.isChecked()
        form_graph_bonds.pushButton_bonds_compare_graphs.setText('Working...')
        form_graph_bonds.pushButton_bonds_compare_graphs.repaint()
        hba.draw_graph(mutations=mutations, compare_to=hba2, draw_edge_occupancy=occupancy, highlight_interregion=interregion, centrality=centrality2, max_centrality=weight, color_dict=color_d, draw_labels=labels, node_factor=0.6)
        form_graph_bonds.close()
        form_graph_bonds.pushButton_bonds_compare_graphs.setText('Compare')
        form_graph_bonds.pushButton_bonds_compare_graphs.repaint()
        
    def save_hbonds_graph():
        color_d={}
        if form.line_bonds_colors.text() != '': color_d = {key.split(':')[0]:key.split(':')[1] for key in form.line_bonds_colors.text().replace(' ','').split(',')}
        labels = form.checkbox_bonds_graph_labels.isChecked()  
        filename = getSaveFileNameWithExt(dialog, 'Save', filter='Portable network graphic (*.png);;Vector graphic (*.eps)')
        if filename == '': return
        form.button_bonds_graph_save.setText('Working...')
        form.button_bonds_graph_save.repaint()
        hba.draw_graph(highlight_interregion=form.checkBox_bonds_interregion.isChecked(), color_dict=color_d, draw_labels=labels, filename=filename)
        form.button_bonds_graph_save.setText('Save')
        form.button_bonds_graph_save.repaint()
        
    def draw_hbonds_occupancy():
        filename = get_filename_hbonds()
        if form_plotting.checkbox_bonds_occupancy_segnames.isChecked():
            draw_hbonds_multi_segname_occupancy()
        else:
            try:
                min_oc = float(form_plotting.line_bonds_occupancy_min.text())
            except:
                min_oc = 0.5
            try:
                max_oc = float(form_plotting.line_bonds_occupancy_max.text())
            except:
                max_oc = 1.0
            try:
                nb_bins = float(form_plotting.line_bonds_occupancy_step.text())
            except:
                nb_bins = 10
            step = (max_oc-min_oc)/nb_bins
            legend = make_legend_list(form_plotting.line_bonds_occupancy_histogram_legend.text())
            f, print_str = hba.draw_occupancy_histogram(min_occupancy=min_oc, max_occupancy=max_oc, occupancy_step=step, legend_text=legend, filename=filename, return_figure=True)
            draw_figure_to_canvas(f, print_str)
        
    def compare_hbonds_occupancy():
        filename = get_filename_hbonds()
        wa2_filename= getOpenFileNameWithExt(dialog_plotting, 'Open', filter='hydrogen bond analysis file (*.hba);;water wire analysis file (*.wwa)')
        if wa2_filename == '': return
        if wa2_filename.split('.')[-1] == 'wa':
            wa2 = WireAnalysis(restore_filename=wa2_filename)
        else:
            wa2 = HbondAnalysis(restore_filename=wa2_filename)
        try:
            min_oc = float(form_plotting.line_bonds_occupancy_min.text())
        except:
            min_oc = 0.5
        try:
            max_oc = float(form_plotting.line_bonds_occupancy_max.text())
        except:
            max_oc = 1.0
        try:
            nb_bins = float(form_plotting.line_bonds_occupancy_step.text())
        except:
            nb_bins = 10
        step = (max_oc-min_oc)/nb_bins
        legend = make_legend_list(form_plotting.line_bonds_occupancy_histogram_legend.text())
        f, print_str = hba.draw_occupancy_histogram(min_occupancy=min_oc, max_occupancy=max_oc, occupancy_step=step, compare_to=wa2, legend_text=legend, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def draw_hbonds_hydration():
        #filename = get_filename_hbonds()
        #try:
        #    cut = float(form_plotting.line_bonds_hydration_min.text())
        #except:
        #    cut = 0.5
        #try:
        #    step = float(form_plotting.line_bonds_hydration_step.text())
        #except:
        #    step = 0.05
        try:
            f, print_str = hba.draw_i4_motif_distribution(return_figure=True)
            draw_figure_to_canvas(f, print_str)
        except:
            error('The analysis has to be initialized with residuewise unchecked')
        
    def draw_hbonds_multi_segname_occupancy():
        filename = get_filename_hbonds()
        colors_temp = str(form_plotting.line_bonds_occupancy_segnames.text()).replace(' ','').split(',')
        colors = {}
        if colors_temp == ['']:
            colors=None
        else:
            try:
                for c in colors_temp:
                    colors[c.split(':')[0]]=c.split(':')[1]
            except:
                colors=None
        try:
            min_oc = float(form_plotting.line_bonds_occupancy_min.text())
        except:
            min_oc = 0.5
        try:
            max_oc = float(form_plotting.line_bonds_occupancy_max.text())
        except:
            max_oc = 1.0
        try:
            nb_bins = float(form_plotting.line_bonds_occupancy_step.text())
        except:
            nb_bins = 10
        step = (max_oc-min_oc)/nb_bins
        f, print_str = hba.draw_multi_segname_occupancy_histogram(min_occupancy=min_oc, max_occupancy=max_oc, occupancy_step=step, colors=colors, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)

    def draw_hbonds_joint_occupancy():
        filename = get_filename_hbonds()
        f, print_str = hba.draw_joint_timeseries(filename=filename, scatter_size=2.0, return_figure=True)
        draw_figure_to_canvas(f, print_str)
    
    def compare_hbonds_joint_occupancy():
        global add_hba
        wa2_filenames= getOpenFileNamesWithExt(dialog_plotting, 'Open', filter='hydrogen bond analysis file (*.hba);;water wire analysis file (*.wwa)')
        if wa2_filenames == '': return
        for wa2_filename in wa2_filenames:
            if wa2_filename.split('.')[-1] == 'wwa':
                wa2 = WireAnalysis(restore_filename=wa2_filename)
            else:
                wa2 = HbondAnalysis(restore_filename=wa2_filename)
            add_hba.append(wa2)
        f, print_str = hba.draw_compare_joint_timeseries(other_paths=add_hba, scatter_size=2.0, return_figure=True)
        add_hba = []
        draw_figure_to_canvas(f, print_str)
        
    def remove_hbonds_joint_occupancy():
        global add_hba
        add_hba = add_hba[:-1]
        f, print_str = hba.draw_compare_joint_timeseries(other_paths=add_hba, scatter_size=2.0, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def draw_hbonds_per_residue():
        filename = get_filename_hbonds()
        #resnames = form_plotting.checkbox_bonds_resnames.isChecked()
        average = form_plotting.checkbox_bonds_per_residue_average.isChecked()
        legend = form_plotting.line_bonds_per_residue_legend.text()
        legend = make_legend_list(legend)
        try:
            ranges = form_plotting.line_bonds_per_residue_range.text()
            ranges = ranges.replace(' ','').split(',')
            ranges_num = []
            for r in ranges:
                if len(r.split('-'))==2:
                    ranges_num += list(range(int(r.split('-')[0]), int(r.split('-')[1])+1))
                elif len(r.split('-'))==1:
                    ranges_num += [int(r.split('-')[0])]
            if ranges_num == []: ranges_num=None
        except:
            ranges_num = None
        f, print_str = hba.draw_connections_per_residue(residues_to_plot=ranges_num, average=average, xtick_resnames=True, legend_text=legend, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
    
    def compare_hbonds_per_residue():
        filename = get_filename_hbonds()
        wa2_filename= getOpenFileNameWithExt(dialog_plotting, 'Open', filter='hydrogen bond analysis file (*.hba);;water wire analysis file (*.wwa)')
        if wa2_filename == '': return
        if wa2_filename.split('.')[-1] == 'wwa':
            wa2 = WireAnalysis(restore_filename=wa2_filename)
        else:
            wa2 = HbondAnalysis(restore_filename=wa2_filename)
        #resnames = form_plotting.checkbox_bonds_resnames.isChecked()
        average = form_plotting.checkbox_bonds_per_residue_average.isChecked()
        legend = form_plotting.line_bonds_per_residue_legend.text()
        legend = make_legend_list(legend)
        if isinstance(legend, list):
            if len(legend) < 2: legend=None
        try:
            ranges = form_plotting.line_bonds_per_residue_range.text().replace(' ','').split(',')
            ranges_num = []
            for r in ranges:
                if len(r.split('-'))==2:
                    ranges_num += list(range(int(r.split('-')[0]), int(r.split('-')[1])+1))
                elif len(r.split('-'))==1:
                    ranges_num += [int(r.split('-')[0])]
            if ranges_num == []: ranges_num = None
        except:
            ranges_num = None
        f, print_str = hba.draw_connections_per_residue(compare_to=wa2, residues_to_plot=ranges_num, average=average, xtick_resnames=True, legend_text=legend, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def draw_hbonds_per_time_interval():
        filename = get_filename_hbonds()
        legend = make_legend_list(form_plotting.line_bonds_time_histogram_legend.text())
        try:
            blocks = int(form_plotting.line_bonds_nb_blocks.text())
        except:
            blocks = 10
        try:
            first = int(form_plotting.line_bonds_time_first.text())
        except:
            first = 0
        try:
            last = int(form_plotting.line_bonds_time_last.text())
        except:
            last = -1
        f, print_str = hba.draw_time_histogram(nb_blocks=blocks+1, first_frame=first, last_frame=last, legend_text=legend, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def compare_hbonds_per_time_interval():
        filename = get_filename_hbonds()
        wa2_filename= getOpenFileNameWithExt(dialog_plotting, 'Open', filter='hydrogen bond analysis file (*.hba);;water wire analysis file (*.wwa)')
        if wa2_filename == '': return
        if wa2_filename.split('.')[-1] == 'wwa':
            wa2 = WireAnalysis(restore_filename=wa2_filename)
        else:
            wa2 = HbondAnalysis(restore_filename=wa2_filename)
        try:
            blocks = int(form_plotting.line_bonds_nb_blocks.text())
        except:
            blocks = 10
        try:
            first = int(form_plotting.line_bonds_time_first.text())
        except:
            first = 0
        try:
            last = int(form_plotting.line_bonds_time_last.text())
        except:
            last = -1
        legend = make_legend_list(form_plotting.line_bonds_time_histogram_legend.text())
        f, print_str = hba.draw_time_histogram(nb_blocks=blocks+1, first_frame=first, last_frame=last, compare_to=wa2, legend_text=legend, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def draw_hbonds_sum_timeseries():
        filename = get_filename_hbonds()
        if form_plotting.checkbox_bonds_timeseries_segnames.isChecked():
            draw_hbonds_multi_segname_timeseries()
        else:
            legend = make_legend_list(form_plotting.line_bonds_sum_legend.text())
            f, print_str = hba.draw_sum_of_connections_timeseries(legend_text=legend, filename=filename, return_figure=True)
            draw_figure_to_canvas(f, print_str)
        
    def draw_hbonds_multi_segname_timeseries():
        filename = get_filename_wires()
        colors_temp = str(form_plotting.line_bonds_timeseries_segnames.text()).replace(' ','').split(',')
        colors = {}
        if colors_temp == ['']:
            colors=None
        else:
            try:
                for c in colors_temp:
                    colors[c.split(':')[0]]=c.split(':')[1]
            except:
                colors=None
        f, print_str = hba.draw_multi_segment_connection_timeseries(colors=colors, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def compare_hbonds_sum_timeseries():
        filename = get_filename_hbonds()
        wa2_filename= getOpenFileNameWithExt(dialog_plotting, 'Open', filter='hydrogen bond analysis file (*.hba);;water wire analysis file (*.wwa)')
        if wa2_filename == '': return
        if wa2_filename.split('.')[-1] == 'wwa':
            wa2 = WireAnalysis(restore_filename=wa2_filename)
        else:
            wa2 = HbondAnalysis(restore_filename=wa2_filename)
        legend = make_legend_list(form_plotting.line_bonds_sum_legend.text())
        if isinstance(legend, list):
            if len(legend) < 2: legend=None
        f, print_str = hba.draw_sum_of_connections_timeseries(compare_to=wa2, legend_text=legend, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def draw_hbonds_heatmap():
        filename = get_filename_hbonds()
        average = form_plotting.checkbox_bonds_heatmap_average.isChecked()
        residuewise = form_plotting.checkbox_bonds_heatmap_residuewise.isChecked()
        ranges = form_plotting.line_bonds_heatmap_ranges.text().split(',')
        names = form_plotting.line_bonds_heatmap_names.text().split(',')
        for i in range(len(names)): names[i]=names[i].strip(' ')
        if residuewise: 
            f, print_str = hba.draw_residue_residue_heatmap(average=average, filename=filename, return_figure=True)
            draw_figure_to_canvas(f, print_str)
            return
        if ranges == ['']: ranges = ['0-'+str(hba.nb_frames)]
        for i in range(len(ranges)):
            r1, r2 = ranges[i].strip(' ').split('-')
            ranges[i] = [int(r1), int(r2)]
        f, print_str = hba.draw_residue_range_heatmap(ranges=ranges, names=names, average=average, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def get_filename_hbonds():
        """
        if form_plotting.group_save.isChecked():
            if str(form_plotting.line_bonds_filename.text()) != '':
                return str(form_plotting.line_bonds_filename.text())
            else: return None
        else: """
        return None
    
    def toggle_hbonds_heatmap_lines():
        residuewise = form_plotting.checkbox_bonds_heatmap_residuewise.isChecked()
        if residuewise:
            form_plotting.line_bonds_heatmap_ranges.setEnabled(False)
            form_plotting.line_bonds_heatmap_names.setEnabled(False)
            form_plotting.label_bonds_ranges.setEnabled(False)
            form_plotting.label_bonds_names.setEnabled(False)
        else:
            form_plotting.line_bonds_heatmap_ranges.setEnabled(True)
            form_plotting.line_bonds_heatmap_names.setEnabled(True)
            form_plotting.label_bonds_ranges.setEnabled(True)
            form_plotting.label_bonds_names.setEnabled(True)
    
    def save_hbonds_plots():
        filename = getSaveFileNameWithExt(dialog_plotting, 'Save', filter='Portable Network Graphics (*.png);;Vector Graphic (*.eps)')
        if not filename: return
        form_plotting.line_bonds_filename.setText(filename)
        
    def save_hbonds_to_file():
        filename = getSaveFileNameWithExt(dialog, 'Save', filter='hydrogen bond analysis file (*.hba)')
        if filename: hba.dump_to_file(filename)
        
    def restore_hbonds_from_file():
        filename = getOpenFileNameWithExt(dialog, 'Open', filter='hydrogen bond analysis file (*.hba)')
        if filename == '': return
        global hba
        form.button_bonds_restore.setText('Working...')
        form.button_bonds_restore.repaint()
        try:
            hba = HbondAnalysis(restore_filename=filename)
            hba._reload_universe()
        except FileNotFoundError:
            estring = 'Couldn`t find the following structure and/or trajectory file(s):\n' 
            if not os.path.isfile(hba._structure): estring += hba._structure + '\n'
            if hba._trajectories is not None: 
                for tra in hba._trajectories: 
                    if not os.path.isfile(tra): estring += tra + '\n'
            error(estring)
            form.button_bonds_restore.setText('Restore State')
            form.button_bonds_restore.repaint()
            return
        form.group_set_bonds.setEnabled(True)
        form.group_bonds_filter.setEnabled(True)
        #form.group_bonds_compute.setEnabled(True)
        form.groupBox_2D_bonds.setEnabled(True)
        form.groupBox_plotting_bonds.setEnabled(True)
        #form.group_bonds_draw.setEnabled(True)
        form.button_bonds_save.setEnabled(True)
        if not hba.residuewise: 
            form_plotting.groupBox.setEnabled(True)
            form.groupBox_filter_backbone.setEnabled(True)
        else: 
            form_plotting.groupBox.setEnabled(False)
            form.groupBox_filter_backbone.setEnabled(False)
            form_plotting.checkbox_bonds_heatmap_residuewise.setEnabled(False)
            form_plotting.checkbox_bonds_heatmap_residuewise.setChecked(False)
            toggle_hbonds_heatmap_lines()
        form.button_bonds_restore.setText('Restore State')
        form.button_bonds_restore.repaint()
        form.line_bonds_structure.setText(hba._structure)
        if hba._trajectories is not None: form.line_bonds_trajectories.setText(', '.join(hba._trajectories))
        else: form.line_bonds_trajectories.setText('')
        if hba._ions_list: form.line_bonds_ions.setText(', '.join(hba._ions_list))
        form.line_bonds_selection.setText(hba._selection)
        form.line_bonds_dist.setText(str(hba.distance))
        form.line_bonds_angle.setText(str(hba.cut_angle))
        start = str(hba._trajectory_slice.start)
        if start == 'None':start = ''
        form.line_bonds_start.setText(start)
        stop = str(hba._trajectory_slice.stop)
        if stop == 'None':stop=''
        form.line_bonds_stop.setText(stop)
        step=str(hba._trajectory_slice.step)
        if step == '1': step=''
        form.line_bonds_step.setText(step)
        add_d = list(hba.donor_names - donor_names_global)
        if add_d: form.line_bonds_donors.setText(', '.join(add_d))
        add_a = list(hba.acceptor_names - acceptor_names_global)
        if add_a: form.line_bonds_acceptors.setText(', '.join(add_a))
        form.line_bonds_missing_residues.setText(str(hba.add_missing_residues))
        form.checkbox_bonds_residuewise.setChecked(hba.residuewise)
        form.checkBox_bonds_donors_without_hydrogen.setChecked(hba._add_donors_without_hydrogen)
        form.checkBox_angle.setChecked(hba.check_angle)
        form.checkBox_bonds_exclude_backbone.setChecked(hba._exclude_backbone_backbone)
        exp, occ, comp, paths, segm, resn, hyd, back = get_applied_filter_hba()
        form.textedit_bonds_active_filters.clear()
        
        if exp is not None: 
            form.textedit_bonds_active_filters.append('Explicit path: '+', '.join(exp))
            form.textedit_bonds_path.append(', '.join(exp))
        if occ is not None: 
            form.textedit_bonds_active_filters.append('Occupancy: '+str(occ))
            form.line_bonds_occupancy.setText(str(occ))
        if comp is not None: 
            form.textedit_bonds_active_filters.append('Connected component: '+comp)
            form.line_bonds_connected_root.setText(comp)
        if paths is not None: 
            form.textedit_bonds_active_filters.append('Paths from '+paths[0]+' to '+paths[1])
            form.line_bonds_path_root.setText(paths[0])
            form.line_bonds_path_goal.setText(paths[1])
        if segm is not None: 
            if segm[1] is not None:
                form.textedit_bonds_active_filters.append('Between segments '+segm[0]+' and '+str(segm[1]))
                form.line_bonds_filter_sega.setText(segm[0])
                form.line_bonds_filter_segb.setText(segm[1])
            else:
                form.textedit_bonds_active_filters.append('Between segment '+segm[0]+' and all')
                form.line_bonds_filter_sega.setText(segm[0])
        if resn is not None: 
            if resn[1] is not None:
                form.textedit_bonds_active_filters.append('Between resnames '+resn[0]+' and '+str(resn[1]))
                form.line_bonds_filter_resa.setText(resn[0])
                form.line_bonds_filter_resb.setText(resn[1])
            else:
                form.textedit_bonds_active_filters.append('Between resname '+resn[0]+' and all')
                form.line_bonds_filter_resa.setText(resn[0])
        if hyd is not None: form.textedit_bonds_active_filters.append('Water in '+str(hyd)+' shells')
        if back is not None: 
            form.textedit_bonds_active_filters.append('Only backbone-backbone H Bonds')
            form.checkBox_bonds_backbone.setChecked(True)
        
    #form_plotting.button_bonds_remove_last.clicked.connect(remove_hbonds_joint_occupancy)
    form_plotting.button_bonds_add_hba.clicked.connect(compare_hbonds_joint_occupancy)
    form_plotting.button_bonds_draw_hydration.clicked.connect(draw_hbonds_hydration)
    form_plotting.button_bonds_draw_occupancy.clicked.connect(draw_hbonds_occupancy)
    form_plotting.button_bonds_compare_occupancy.clicked.connect(compare_hbonds_occupancy)
    form_plotting.button_bonds_draw_bonds_per_residue.clicked.connect(draw_hbonds_per_residue)
    form_plotting.button_bonds_compare_bonds_per_residue.clicked.connect(compare_hbonds_per_residue)
    form_plotting.button_bonds_draw_joint_barcode.clicked.connect(draw_hbonds_joint_occupancy)
    form_plotting.button_bonds_draw_bonds_timeseries.clicked.connect(draw_hbonds_sum_timeseries)
    form_plotting.button_bonds_compare_bonds_timeseries.clicked.connect(compare_hbonds_sum_timeseries)
    form_plotting.button_bonds_draw_bonds_per_time_interval.clicked.connect(draw_hbonds_per_time_interval)
    form_plotting.button_bonds_compare_bonds_per_time_interval.clicked.connect(compare_hbonds_per_time_interval)
    form_plotting.button_bonds_draw_heatmap.clicked.connect(draw_hbonds_heatmap)
    form_plotting.checkbox_bonds_heatmap_residuewise.stateChanged.connect(toggle_hbonds_heatmap_lines)
    form_plotting.checkbox_bonds_occupancy_segnames.stateChanged.connect(toggle_bonds_occupancy_segments)
    form_plotting.checkbox_bonds_timeseries_segnames.stateChanged.connect(toggle_bonds_timeseries_segments)
    form.button_bonds_open_plotting.clicked.connect(load_plotting_bonds)
    form.button_bonds_draw_graph.clicked.connect(draw_bonds_graph)
    #form.button_bonds_compute_residence.clicked.connect(compute_bonds_residence)
    #form.button_bonds_compute_i4.clicked.connect(compute_i4)
    #form.button_bonds_compute_occupancy.clicked.connect(compute_bonds_joint)
    form.button_bonds_filter.clicked.connect(filter_bonds)
    form.button_bonds_clear.clicked.connect(clear_all_filters_bonds)
    form.button_bonds_structure.clicked.connect(browse_filename_structure)
    form.button_bonds_trajectories.clicked.connect(browse_filename_trajectories)
    form.button_bonds_init.clicked.connect(init_bonds_analysis)
    form.bonds_set.clicked.connect(set_bonds)
    form.button_bonds_save.clicked.connect(save_hbonds_to_file)
    form.button_bonds_restore.clicked.connect(restore_hbonds_from_file)
    #form.button_bonds_graph_save.clicked.connect(save_hbonds_graph)
    form.checkBox_angle.stateChanged.connect(toggle_bonds_angle)
    form.pushButton_bonds_save_data.clicked.connect(save_bonds_data)
    form.pushButton_bonds_advanced.clicked.connect(load_advanced_bonds)
    form_graph_bonds.pushButton_bonds_save_data.clicked.connect(save_bonds_data_advanced)
    form_graph_bonds.pushButton_draw_graph.clicked.connect(draw_bonds_graph_advanced)
    form_graph_bonds.pushButton_bonds_compare_graphs.clicked.connect(draw_bonds_graph_compare)
    
    def browse_wire_structure():
        filename = getOpenFileNameWithExt(
            dialog, 'Open', filter='PSF File (*.psf);;PDB File (*.pdb);;All Files (*.*)')
        if filename:
            form.line_wire_structure.setText(filename)
    
    def browse_wire_trajectories():
        filename = getOpenFileNamesWithExt(
            dialog, 'Open', filter='DCD File (*.dcd);;All Files (*.*)')
        if filename:
            fnames = ', '.join(filename)
            form.line_wire_trajectories.setText(fnames)
    
    def init_wire_analysis():
        global wa
        form.button_wire_init.setText('Working...')
        form.button_wire_init.repaint()
        ad = [donor for donor in form.line_wire_donors.text().replace(' ','').split(',')]
        aa = [acceptor for acceptor in form.line_wire_acceptors.text().replace(' ','').split(',')]
        form.group_bonds_filter_2.setEnabled(False)
        form.group_bonds_filter_2.repaint()
        #form.group_wire_compute.setEnabled(True)
        form.groupBox_plotting_wire.setEnabled(False)
        form.groupBox_2D_wire.setEnabled(False)
        form.groupBox_2D_wire.repaint()
        #form.group_wire_draw.setEnabled(True)
        form.button_wire_save.setEnabled(False)
        form.button_wire_save.repaint()
        form.group_wire_set.setEnabled(False)
        form.group_wire_set.repaint()
        form.textedit_bonds_active_filters_2.clear()
        form.repaint()
        set_params_to_default_wire()
        try:
            step = form.line_wire_step.text()
            stop = form.line_wire_stop.text()
            start = form.line_wire_start.text()
            if step == '': step=1
            else: step = int(step)
            if stop == '': stop=None
            else: stop = int(stop)
            if start == '': start=None
            else: start = int(start)
            trajectory = form.line_wire_trajectories.text().split(', ')
            if trajectory == ['']: trajectory = None
            donors_without = form.checkBox_wires_donors_without_hydrogen.isChecked()
            try:
                wa = WireAnalysis(form.line_wire_selection.text(), form.line_wire_structure.text(),
                                    trajectory, float(form.line_wire_dist.text()),
                                    float(form.line_wire_angle.text()), start=start,
                                    stop=stop, step=step, residuewise=form.checkbox_wires_residuewise.isChecked(),
                                    additional_donors=ad, additional_acceptors=aa, add_donors_without_hydrogen=donors_without)
                form.group_wire_set.setEnabled(True)
                try:wa.add_missing_residues = int(form.line_wire_missing_residues.text())
                except: error('Missing residues number not set. Cant convert input to integer.')
            except AssertionError:
                error('Couldnt initialize analysis. Is water present?')
                form.button_wire_init.setText('Initialize')
                form.button_wire_init.repaint()
                return
            if form.checkbox_wires_residuewise.isChecked():
                form_plotting.checkbox_wires_heatmap_residuewise.setEnabled(False)
                form_plotting.checkbox_wires_heatmap_residuewise.setChecked(False)
                toggle_wire_heatmap_lines()
            form.button_wire_init.setText('Initialize')
            form.button_wire_init.repaint()
        except: 
            error('One or more input fields in the Initialize area are of unexpected type or empty.')
            form.button_wire_init.setText('Initialize')
            form.button_wire_init.repaint()
        
    def set_wires():
        try:
            step = form.line_wire_step.text()
            stop = form.line_wire_stop.text()
            start = form.line_wire_start.text()
            if step == '': step=1
            else: step = int(step)
            if stop == '': stop=None
            else: stop = int(stop)
            if start == '': start=None
            else: start = int(start)
            dist = float(form.line_wire_dist.text())
            angle = float(form.line_wire_angle.text())
            check_angle = form.checkBox_wires_angle.isChecked()
            direct_bonds = form.checkBox_wires_allow_direct_bonds.isChecked()
            if form.checkBox_wires_donors_without_hydrogen.isChecked() and check_angle: 
                error('Donors without hydrogen can only be used with no angle check! Proceeding with angle check set to false')
                form.checkBox_wires_angle.setChecked(False)
                check_angle = False
            wa.distance = dist
            wa.cut_angle = angle
            wa.check_angle = check_angle
            wa._trajectory_slice = slice(start if isinstance(start, int) else None, stop if isinstance(stop, int) else None, step)
            wa.nb_frames = len([0 for i in wa._universe.trajectory[wa._trajectory_slice]])
        except:
            error('One or more input fields in the Search area are of unexpected type.')
        try:
            max_water = int(form.line_wire_max_water.text())
        except:
            max_water = 5 
            form.line_wire_max_water.setText('5')
        form.button_wire_set.setText('Working...')
        form.button_wire_set.repaint()
        if form.radio_wire_water.isChecked():
            wa.set_water_wires_csr(max_water, allow_direct_bonds=direct_bonds)
        else:
            wa.set_water_wires(max_water, allow_direct_bonds=direct_bonds)
        form.button_wire_set.setText('Search')
        clear_all_filters_wires()
        form.group_bonds_filter_2.setEnabled(True)
        #form.group_wire_compute.setEnabled(True)
        form.groupBox_plotting_wire.setEnabled(True)
        form.groupBox_2D_wire.setEnabled(True)
        #form.group_wire_draw.setEnabled(True)
        form.button_wire_save.setEnabled(True)
    
    def clear_all_filters_wires():
        form.textedit_bonds_active_filters_2.clear()
        form.line_bonds_occupancy_2.setText('')
        form.textedit_bonds_path_2.setText('')
        form.line_bonds_connected_root_2.setText('')
        form.line_bonds_path_root_2.setText('')
        form.line_bonds_path_goal_2.setText('')
        form.line_bonds_filter_sega_2.setText('')
        form.line_bonds_filter_segb_2.setText('')
        form.line_bonds_filter_resa_2.setText('')
        form.line_bonds_filter_resb_2.setText('')
        form.line_wires_order.setText('')
        filter_wires()
        form.repaint()
        
    def set_params_to_default_wire():
        form.line_wire_start.clear()
        form.line_wire_stop.clear()
        form.line_wire_step.clear()
        form.line_wire_dist.setText('3.5')
        form.line_wire_angle.setText('60')
        form.repaint()
    
    def filter_wires():
        occupancy_cut = form.line_bonds_occupancy_2.text()
        connected_root = form.line_bonds_connected_root_2.text()
        path_root = form.line_bonds_path_root_2.text()
        path_goal = form.line_bonds_path_goal_2.text()
        path_explicit = form.textedit_bonds_path_2.toPlainText()
        sega = form.line_bonds_filter_sega_2.text()
        segb = form.line_bonds_filter_segb_2.text()
        resa = form.line_bonds_filter_resa_2.text()
        resb = form.line_bonds_filter_resb_2.text()
        form.textedit_bonds_active_filters_2.clear()
        use_filtered = False
        try:
            order = [int(a) for a in form.line_wires_order.text().replace(' ','').split(',') if a!='']
            if len(order) == 0: order=[3,1,2,4,5] 
        except:
            error('order must be integers separated by commas or empty')
            return
        for a in order:
            if path_explicit != '' and a==3:
                try:
                    args = [node for node in path_explicit.replace(' ','').strip('\n').split(',')]
                    wa.filter_single_path(args, use_filtered)
                    form.textedit_bonds_active_filters_2.append('Explicit path: '+', '.join(args))
                    use_filtered = True
                except: error('Either one or more of the nodes specified are not in the graph, or the nodes dont build a path in that order!')
            if occupancy_cut != '' and a==1:
                try: 
                    occupancy_cut = float(occupancy_cut)
                    wa.filter_occupancy(occupancy_cut, use_filtered)
                    form.textedit_bonds_active_filters_2.append('Occupancy: '+str(occupancy_cut))
                    use_filtered = True
                except: error('min. value has to be a float!')
            if connected_root != '' and a==2:
                try:
                    wa.filter_connected_component(connected_root, use_filtered)
                    form.textedit_bonds_active_filters_2.append('Connected component: '+connected_root)
                    use_filtered = True
                except: error('The specified root is not in the graph!')
            if (path_root != '' and path_goal != '') and a==4:
                try:
                    if form.radioButton_wires_number_wires.isChecked():
                        wa.filter_all_paths(path_root, path_goal, use_filtered)
                        form.textedit_bonds_active_filters_2.append('Shortest paths from '+path_root+' to '+path_goal)
                        use_filtered = True
                    else:
                        wa.filter_minmal_bonds_path(path_root, path_goal, use_filtered)
                        form.textedit_bonds_active_filters_2.append('Minimal bonds paths from '+path_root+' to '+path_goal)
                        use_filtered = True
                except: error('root and/or target are not in the graph!')
            if (sega != '' or resa != '') and a==5:
                try:
                    if sega != '':
                        if segb == '': segb=None
                        wa.filter_between_segnames(sega, segb, use_filtered)
                        if segb is not None: form.textedit_bonds_active_filters_2.append('Between segments '+sega+' and '+segb)
                        else: form.textedit_bonds_active_filters_2.append('Between segment '+sega+' and all')
                        use_filtered = True
                    if resa != '':
                        if resb == '': resb=None
                        wa.filter_between_resnames(resa, resb, use_filtered)
                        if resb is not None: form.textedit_bonds_active_filters_2.append('Between resnames '+resa+' and '+resb)
                        else: form.textedit_bonds_active_filters_2.append('Between resname '+resa+' and all')
                        use_filtered = True
                except:
                    pass
                
        if len(wa.filtered_results) == 0:
            error('This combination of filters leads to an empty set of results!')
            use_filtered = False
        if not use_filtered: 
            wa.filtered_graph = wa.initial_graph
            wa.filtered_results = wa.initial_results
            form.textedit_bonds_active_filters_2.clear()
            
    def compute_wires_joint():
        joint_occupancy = wa.compute_joint_occupancy()
        form.textedit_wire_results.clear()
        if joint_occupancy.mean()>0.0:
            form.textedit_wire_results.append('Joint occupancy: '+str(np.round(joint_occupancy.mean(), 2))+', first frame: '+str(np.argmax(joint_occupancy)))
        else:
            form.textedit_wire_results.append('Joint occupancy: 0.0')
            
    def compute_average():
        import numpy as _np
        projection = wa.compute_average_water_per_wire()
        projection_printable = ''
        for key in projection:
            projection_printable += key+ '.'*(45-len(key)) + ' | '+ str(_np.round(projection[key],2))+'\n'
        form.textedit_wire_results.clear()
        form.textedit_wire_results.setText(projection_printable)
        
    def compute_wire_projection():
        import numpy as _np
        projection = wa.compute_wire_projection()
        projection_printable = ''
        for key in projection:
            projection_printable += key+ '.'*(45-len(key)) + ' | '+ str(_np.round(_np.mean(projection[key][projection[key]!=-1]),2))+'\n'
        form.textedit_wire_results.clear()
        form.textedit_wire_results.setText(projection_printable)
    
    def draw_wire_graph():
        form.button_wire_draw_graph.setText('Working...')
        form.button_wire_draw_graph.repaint()
        wa.draw_graph(draw_edge_occupancy=False, draw_labels=True, node_factor=0.6)
        form.button_wire_draw_graph.setText('Draw')
        form.button_wire_draw_graph.repaint()
        
    def draw_wire_graph_advanced():
        color_d={}
        if form_graph_wires.line_bonds_colors.text() != '': color_d = {key.split(':')[0]:key.split(':')[1] for key in form_graph_wires.line_bonds_colors.text().replace(' ','').split(',')}
        centrality = None
        weight = None
        if form_graph_wires.lineEdit_bonds_color_weight.text() != '': 
            try: weight = float(form_graph_wires.lineEdit_bonds_color_weight.text())
            except: pass
        if form_graph_wires.radioButton_bonds_betweenness.isChecked():
            centrality = wa.compute_centrality(centrality_type='betweenness', weight=None)
        elif form_graph_wires.radioButton_bonds_degree.isChecked():
            centrality = wa.compute_centrality(centrality_type='degree', normalize=True, weight=None)
        labels = form_graph_wires.checkBox_bonds_graph_labels.isChecked()
        interregion = form_graph_wires.checkBox_bonds_interregion.isChecked()
        occupancy = form_graph_wires.checkBox_bonds_occupancy.isChecked()
        form_graph_wires.pushButton_draw_graph.setText('Working...')
        form_graph_wires.pushButton_draw_graph.repaint()
        wa.draw_graph(draw_edge_occupancy=occupancy, highlight_interregion=interregion, centrality=centrality, max_centrality=weight, color_dict=color_d, draw_labels=labels, node_factor=0.6)
        form_graph_wires.close()
        form_graph_wires.pushButton_draw_graph.setText('Draw')
        form_graph_wires.pushButton_draw_graph.repaint()
        
    def draw_wire_graph_compare():
        filename = getOpenFileNameWithExt(dialog_graph_bonds, 'Open', filter='water wire analysis file (*.wwa);;hydrogen bond analysis file (*.hba)')
        if not filename: return
        hba2 = HbondAnalysis(restore_filename=filename)
        mutations=None
        if form_graph_wires.lineEdit_mutations.text() != '':
            mutations = form_graph_wires.lineEdit_mutations.text().replace(' ', '').split(',')
        color_d={}
        if form_graph_wires.line_bonds_colors.text() != '': color_d = {key.split(':')[0]:key.split(':')[1] for key in form_graph_wires.line_bonds_colors.text().replace(' ','').split(',')}
        centrality = None
        weight = None
        if form_graph_wires.lineEdit_bonds_color_weight.text() != '': 
            try: weight = float(form_graph_wires.lineEdit_bonds_color_weight.text())
            except: pass
        if form_graph_wires.radioButton_bonds_betweenness.isChecked():
            if form_graph_wires.checkBox_bonds_weight_centrality.isChecked(): centrality = wa.compute_centrality(centrality_type='betweenness', weight=None)
            centrality2 = wa.compute_centrality(centrality_type='betweenness', weight=None)
        elif form_graph_wires.radioButton_bonds_degree.isChecked():
            if form_graph_wires.checkBox_bonds_weight_centrality.isChecked(): centrality = wa.compute_centrality(centrality_type='degree', normalize=True, weight=None)
            centrality2 = wa.compute_centrality(centrality_type='degree', weight=None)
        if form_graph_wires.checkBox_bonds_weight_centrality.isChecked() and (form_graph_bonds.radioButton_bonds_betweenness.isChecked() or form_graph_bonds.radioButton_bonds_degree.isChecked()):
            ma = np.max([np.max(list(centrality2.values())), np.max(list(centrality.values()))])
            for node in centrality2:
                centrality2[node] /= ma
        labels = form_graph_wires.checkBox_bonds_graph_labels.isChecked()
        interregion = form_graph_wires.checkBox_bonds_interregion.isChecked()
        occupancy = form_graph_wires.checkBox_bonds_occupancy.isChecked()
        form_graph_wires.pushButton_bonds_compare_graphs.setText('Working...')
        form_graph_wires.pushButton_bonds_compare_graphs.repaint()
        wa.draw_graph(mutations=mutations, compare_to=hba2, draw_edge_occupancy=occupancy, highlight_interregion=interregion, centrality=centrality, max_centrality=weight, color_dict=color_d, draw_labels=labels, node_factor=0.6)
        form_graph_wires.close()
        form_graph_wires.pushButton_bonds_compare_graphs.setText('Compare')
        form_graph_wires.pushButton_bonds_compare_graphs.repaint()
    
    def save_wire_graph():
        color_d={}
        if form.line_wire_colors.text() != '': color_d = {key.split(':')[0]:key.split(':')[1] for key in form.line_wire_colors.text().replace(' ','').split(',')}
        labels = form.checkbox_wire_graph_labels.isChecked()  
        filename = getSaveFileNameWithExt(dialog, 'Save', filter='Portable network graphic (*.png);;Vector Graphic (*.eps)')
        if filename == '': return
        form.button_wire_save_graph.setText('Working...')
        form.button_wire_save_graph.repaint()
        wa.draw_graph(highlight_interregion=form.checkBox_wire_interregion.isChecked(), color_dict=color_d, draw_labels=labels, filename=filename)
        form.button_wire_save_graph.setText('Save')
        form.button_wire_save_graph.repaint()
    
    def toggle_wire_occupancy_segments():
        form_plotting.line_wires_occupancy_segnames.setEnabled(form_plotting.checkbox_wires_occupancy_segnames.isChecked())
        form_plotting.line_wires_occupancy_histogram_legend.setEnabled(not form_plotting.checkbox_wires_occupancy_segnames.isChecked())
        form_plotting.label_18.setEnabled(not form_plotting.checkbox_wires_occupancy_segnames.isChecked())
        form_plotting.button_wires_compare_occupancy.setEnabled(not form_plotting.checkbox_wires_occupancy_segnames.isChecked())
    
    def toggle_wire_timeseries_segments():
        form_plotting.line_wires_timeseries_segnames.setEnabled(form_plotting.checkbox_wires_timeseries_segnames.isChecked())
        form_plotting.line_wires_sum_legend.setEnabled(not form_plotting.checkbox_wires_timeseries_segnames.isChecked())
        form_plotting.label_7.setEnabled(not form_plotting.checkbox_wires_timeseries_segnames.isChecked())
        form_plotting.button_wires_compare_bonds_timeseries.setEnabled(not form_plotting.checkbox_wires_timeseries_segnames.isChecked())
    
    def draw_wire_occupancy():
        filename = get_filename_wires()
        if form_plotting.checkbox_wires_occupancy_segnames.isChecked(): 
            draw_wire_multi_segname_occupancy()
        else:
            try:
                min_oc = float(form_plotting.line_wires_occupancy_min.text())
            except:
                min_oc = 0.5
            try:
                max_oc = float(form_plotting.line_wires_occupancy_max.text())
            except:
                max_oc = 1.0
            try:
                nb_bins = float(form_plotting.line_wires_occupancy_step.text())
            except:
                nb_bins = 10
            step = (max_oc-min_oc)/nb_bins
            legend = make_legend_list(form_plotting.line_wires_occupancy_histogram_legend.text())
            f, print_str = wa.draw_occupancy_histogram(min_occupancy=min_oc, max_occupancy=max_oc, occupancy_step=step, legend_text=legend, filename=filename, return_figure=True)
            draw_figure_to_canvas(f, print_str)
        
    def compare_wire_occupancy():
        filename = get_filename_wires()
        wa2_filename= getOpenFileNameWithExt(dialog_plotting, 'Open', filter='water wire analysis file (*.wwa);;hydrogen bond analysis file (*.hba)')
        if wa2_filename == '': return
        if wa2_filename.split('.')[-1] == 'wwa':
            wa2 = WireAnalysis(restore_filename=wa2_filename)
        else:
            wa2 = HbondAnalysis(restore_filename=wa2_filename)
        try:
            min_oc = float(form_plotting.line_wires_occupancy_min.text())
        except:
            min_oc = 0.5
        try:
            max_oc = float(form_plotting.line_wires_occupancy_max.text())
        except:
            max_oc = 1.0
        try:
            nb_bins = float(form_plotting.line_wires_occupancy_step.text())
        except:
            nb_bins = 10
        step = (max_oc-min_oc)/nb_bins
        legend = make_legend_list(form_plotting.line_wires_occupancy_histogram_legend.text())
        f, print_str = wa.draw_occupancy_histogram(min_occupancy=min_oc, max_occupancy=max_oc, occupancy_step=step, compare_to=wa2, legend_text=legend, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def draw_wire_multi_segname_occupancy():
        filename = get_filename_wires()
        colors_temp = str(form_plotting.line_wires_occupancy_segnames.text()).replace(' ','').split(',')
        colors = {}
        if colors_temp == ['']:
            colors=None
        else:
            try:
                for c in colors_temp:
                    colors[c.split(':')[0]]=c.split(':')[1]
            except:
                colors=None
        try:
            min_oc = float(form_plotting.line_wires_occupancy_min.text())
        except:
            min_oc = 0.5
        try:
            max_oc = float(form_plotting.line_wires_occupancy_max.text())
        except:
            max_oc = 1.0
        try:
            nb_bins = float(form_plotting.line_wires_occupancy_step.text())
        except:
            nb_bins = 10
        step = (max_oc-min_oc)/nb_bins
        f, print_str = wa.draw_multi_segname_occupancy_histogram(min_occupancy=min_oc, max_occupancy=max_oc, occupancy_step=step, colors=colors, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def draw_wire_joint_occupancy():
        filename = get_filename_wires()
        f, print_str = wa.draw_joint_timeseries(filename=filename, scatter_size=2.0, return_figure=True)
        draw_figure_to_canvas(f, print_str)
    
    def draw_wire_water_in_wire():
        resa = form_plotting.line_wires_in_wire_res_a.text()
        resb = form_plotting.line_wires_in_wire_res_b.text()
        if resa != '' and resb != '':
            try:
                f, print_str = wa.draw_water_timeseries(resa=resa , resb=resb , scatter_size=2.0, return_figure=True)
                draw_figure_to_canvas(f, print_str)
            except KeyError:
                error('No connection between the specified residues found')
        else:
            error('Specify two connected residues!')
    
    def compare_wire_joint_occupancy():
        global add_wa
        wa2_filenames= getOpenFileNamesWithExt(dialog_plotting, 'Open', filter='water wire analysis file (*.wwa);;hydrogen bond analysis file (*.hba)')
        if wa2_filenames == '': return
        for wa2_filename in wa2_filenames:
            if wa2_filename.split('.')[-1] == 'wwa':
                wa2 = WireAnalysis(restore_filename=wa2_filename)
            else:
                wa2 = HbondAnalysis(restore_filename=wa2_filename)
            add_wa.append(wa2)
        f, print_str = wa.draw_compare_joint_timeseries(other_paths=add_wa, scatter_size=2.0, return_figure=True)
        add_wa=[]
        draw_figure_to_canvas(f, print_str)
        
    def remove_wire_joint_occupancy():
        global add_wa
        add_wa = add_wa[:-1]
        f, print_str = wa.draw_compare_joint_timeseries(other_paths=add_wa, return_figure=True)
        draw_figure_to_canvas(f, print_str)
    
    def draw_wire_per_residue():
        filename = get_filename_wires()
        #resnames = form_plotting.checkbox_bonds_resnames.isChecked()
        average = form_plotting.checkbox_wires_per_residue_average.isChecked()
        legend = form_plotting.line_wires_per_residue_legend.text()
        legend = make_legend_list(legend)
        try:
            ranges = form_plotting.line_wires_per_residue_range.text()
            ranges = ranges.replace(' ','').split(',')
            ranges_num = []
            for r in ranges:
                if len(r.split('-'))==2:
                    ranges_num += list(range(int(r.split('-')[0]), int(r.split('-')[1])+1))
                elif len(r.split('-'))==1:
                    ranges_num += [int(r.split('-')[0])]
            if ranges_num == []: ranges_num=None
        except:
            ranges_num = None
        f, print_str = wa.draw_connections_per_residue(residues_to_plot=ranges_num, average=average, xtick_resnames=True, legend_text=legend, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
    
    def compare_wire_per_residue():
        filename = get_filename_wires()
        wa2_filename= getOpenFileNameWithExt(dialog_plotting, 'Open', filter='water wire analysis file (*.wwa);;hydrogen bond analysis file (*.hba)')
        if wa2_filename == '': return
        if wa2_filename.split('.')[-1] == 'wwa':
            wa2 = WireAnalysis(restore_filename=wa2_filename)
        else:
            wa2 = HbondAnalysis(restore_filename=wa2_filename)
        #resnames = form_plotting.checkbox_bonds_resnames.isChecked()
        average = form_plotting.checkbox_wires_per_residue_average.isChecked()
        legend = form_plotting.line_wires_per_residue_legend.text()
        legend = make_legend_list(legend)
        if isinstance(legend, list):
            if len(legend) < 2: legend=None
        try:
            ranges = form_plotting.line_wires_per_residue_range.text().replace(' ','').split(',')
            ranges_num = []
            for r in ranges:
                if len(r.split('-'))==2:
                    ranges_num += list(range(int(r.split('-')[0]), int(r.split('-')[1])+1))
                elif len(r.split('-'))==1:
                    ranges_num += [int(r.split('-')[0])]
            if ranges_num == []: ranges_num = None
        except:
            ranges_num = None
        f, print_str = wa.draw_connections_per_residue(compare_to=wa2, residues_to_plot=ranges_num, average=average, xtick_resnames=True, legend_text=legend, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def draw_wire_per_time_interval():
        filename = get_filename_wires()
        legend = make_legend_list(form_plotting.line_wires_time_histogram_legend.text())
        try:
            blocks = int(form_plotting.line_wires_nb_blocks.text())
        except:
            blocks = 10
        try:
            first = int(form_plotting.line_wires_time_first.text())
        except:
            first = 0
        try:
            last = int(form_plotting.line_wires_time_last.text())
        except:
            last = -1
        f, print_str = wa.draw_time_histogram(nb_blocks=blocks+1, first_frame=first, last_frame=last, legend_text=legend, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def compare_wire_per_time_interval():
        filename = get_filename_wires()
        wa2_filename= getOpenFileNameWithExt(dialog_plotting, 'Open', filter='water wire analysis file (*.wwa);;hydrogen bond analysis file (*.hba)')
        if wa2_filename == '': return
        if wa2_filename.split('.')[-1] == 'wwa':
            wa2 = WireAnalysis(restore_filename=wa2_filename)
        else:
            wa2 = HbondAnalysis(restore_filename=wa2_filename)
        try:
            blocks = int(form_plotting.line_wires_nb_blocks.text())
        except:
            blocks = 10
        try:
            first = int(form_plotting.line_wires_time_first.text())
        except:
            first = 0
        try:
            last = int(form_plotting.line_wires_time_last.text())
        except:
            last = -1
        legend = make_legend_list(form_plotting.line_wires_time_histogram_legend.text())
        f, print_str = wa.draw_time_histogram(nb_blocks=blocks+1, first_frame=first, last_frame=last, compare_to=wa2, legend_text=legend, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def draw_wire_sum_timeseries():
        filename = get_filename_wires()
        if form_plotting.checkbox_wires_timeseries_segnames.isChecked():
            draw_wire_multi_segname_timeseries()
        else:
            legend = make_legend_list(form_plotting.line_wires_sum_legend.text())
            f, print_str = wa.draw_sum_of_connections_timeseries(legend_text=legend, filename=filename, return_figure=True)
            draw_figure_to_canvas(f, print_str)
        
    def draw_wire_multi_segname_timeseries():
        filename = get_filename_wires()
        colors_temp = str(form_plotting.line_wires_timeseries_segnames.text()).replace(' ','').split(',')
        colors = {}
        if colors_temp == ['']:
            colors=None
        else:
            try:
                for c in colors_temp:
                    colors[c.split(':')[0]]=c.split(':')[1]
            except:
                colors=None
        f, print_str = wa.draw_multi_segment_connection_timeseries(colors=colors, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def compare_wire_sum_timeseries():
        filename = get_filename_wires()
        wa2_filename= getOpenFileNameWithExt(dialog_plotting, 'Open', filter='water wire analysis file (*.wwa);;hydrogen bond analysis file (*.hba)')
        if wa2_filename == '': return
        if wa2_filename.split('.')[-1] == 'wwa':
            wa2 = WireAnalysis(restore_filename=wa2_filename)
        else:
            wa2 = HbondAnalysis(restore_filename=wa2_filename)
        legend = make_legend_list(form_plotting.line_wires_sum_legend.text())
        if isinstance(legend, list):
            if len(legend) < 2: legend=None
        f, print_str = wa.draw_sum_of_connections_timeseries(compare_to=wa2, legend_text=legend, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def draw_wire_heatmap():
        filename = get_filename_wires()
        average = form_plotting.checkbox_wires_heatmap_average.isChecked()
        residuewise = form_plotting.checkbox_wires_heatmap_residuewise.isChecked()
        ranges = form_plotting.line_wires_heatmap_ranges.text().split(',')
        names = form_plotting.line_wires_heatmap_names.text().split(',')
        for i in range(len(names)): names[i]=names[i].strip(' ')
        if residuewise: 
            f, print_str = wa.draw_residue_residue_heatmap(average=average, filename=filename, return_figure=True)
            draw_figure_to_canvas(f, print_str)
            return
        if ranges == ['']: ranges = ['0-'+str(wa.nb_frames)]
        for i in range(len(ranges)):
            r1, r2 = ranges[i].strip(' ').split('-')
            ranges[i] = [int(r1), int(r2)]
        f, print_str = wa.draw_residue_range_heatmap(ranges=ranges, names=names, average=average, filename=filename, return_figure=True)
        draw_figure_to_canvas(f, print_str)
        
    def toggle_wire_heatmap_lines():
        residuewise = form_plotting.checkbox_wires_heatmap_residuewise.isChecked()
        if residuewise:
            form_plotting.line_wires_heatmap_ranges.setEnabled(False)
            form_plotting.line_wires_heatmap_names.setEnabled(False)
            form_plotting.label_wires_names.setEnabled(False)
            form_plotting.label_wires_ranges.setEnabled(False)
        else:
            form_plotting.line_wires_heatmap_ranges.setEnabled(True)
            form_plotting.line_wires_heatmap_names.setEnabled(True)
            form_plotting.label_wires_names.setEnabled(True)
            form_plotting.label_wires_ranges.setEnabled(True)
    
    def get_filename_wires():
        """
        if form_plotting.group_save.isChecked():
            if str(form_plotting.line_bonds_filename.text()) != '':
                return str(form_plotting.line_bonds_filename.text())
            else: return None
        else: """
        return None
    
    def save_wires_plots():
        filename = getSaveFileNameWithExt(dialog_plotting, 'Save', filter='Portable Network Graphics (*.png);;Vector Graphic (*.eps)')
        if not filename: return
        form_plotting.line_wires_filename.setText(filename)
    
    def save_wires_to_file():
        filename = getSaveFileNameWithExt(dialog, 'Save', filter='water wire analysis file (*.wwa)')
        if filename: wa.dump_to_file(filename)
        
    def restore_wires_from_file():
        filename = getOpenFileNameWithExt(dialog, 'Open', filter='water wire analysis file (*.wwa)')
        if filename == '': return
        global wa
        form.button_wire_restore.setText('Working...')
        form.button_wire_restore.repaint()
        try:
            wa = WireAnalysis(restore_filename=filename)
            wa._reload_universe()
        except FileNotFoundError:
            estring = 'Couldn`t find the following structure and/or trajectory file(s):\n' 
            if not os.path.isfile(wa._structure): estring += wa._structure + '\n'
            if wa._trajectories is not None: 
                for tra in wa._trajectories: 
                    if not os.path.isfile(tra): estring += tra + '\n'
            error(estring)
            form.button_wire_restore.setText('Restore State')
            form.button_wire_restore.repaint()
            return
        form.group_wire_set.setEnabled(True)
        form.group_bonds_filter_2.setEnabled(True)
        #form.group_wire_compute.setEnabled(True)
        form.groupBox_2D_wire.setEnabled(True)
        form.groupBox_plotting_wire.setEnabled(True)
        #form.group_wire_draw.setEnabled(True)
        form.button_wire_save.setEnabled(True)
        form.button_wire_restore.setText('Restore State')
        form.button_wire_restore.repaint()
        form.line_wire_structure.setText(wa._structure)
        if wa.residuewise:
            form_plotting.checkbox_wires_heatmap_residuewise.setEnabled(False)
            form_plotting.checkbox_wires_heatmap_residuewise.setChecked(False)
            toggle_wire_heatmap_lines()
        else:
            form_plotting.checkbox_wires_heatmap_residuewise.setEnabled(True)
        if wa._trajectories is not None: form.line_wire_trajectories.setText(', '.join(wa._trajectories))
        else: form.line_wire_trajectories.setText('')
        form.line_wire_selection.setText(wa._selection)
        form.line_wire_dist.setText(str(wa.distance))
        form.line_wire_angle.setText(str(wa.cut_angle))
        start = str(wa._trajectory_slice.start)
        if start == 'None':start = ''
        form.line_wire_start.setText(start)
        stop = str(wa._trajectory_slice.stop)
        if stop == 'None':stop=''
        form.line_wire_stop.setText(stop)
        step=str(wa._trajectory_slice.step)
        if step == '1': step=''
        form.line_wire_step.setText(step)
        add_d = list(wa.donor_names - donor_names_global)
        if add_d: form.line_wire_donors.setText(', '.join(add_d))
        add_a = list(wa.acceptor_names - acceptor_names_global)
        if add_a: form.line_wire_acceptors.setText(', '.join(add_a))
        form.checkbox_wires_residuewise.setChecked(wa.residuewise)
        form.checkBox_wires_donors_without_hydrogen.setChecked(wa._add_donors_without_hydrogen)
        form.checkBox_wires_angle.setChecked(wa.check_angle)
        form.checkBox_wires_allow_direct_bonds.setChecked(wa._allow_direct_bonds)
        exp, occ, comp, paths, segm, resn, avg_least = get_applied_filter_wa()
        form.textedit_bonds_active_filters_2.clear()
        form.line_wire_missing_residues.setText(str(wa.add_missing_residues))
        if exp is not None: 
            form.textedit_bonds_active_filters_2.append('Explicit path: '+', '.join(exp))
            form.textedit_bonds_path_2.append(', '.join(exp))
        if occ is not None: 
            form.textedit_bonds_active_filters_2.append('Occupancy: '+str(occ))
            form.line_bonds_occupancy_2.setText(str(occ))
        if comp is not None: 
            form.textedit_bonds_active_filters_2.append('Connected component: '+comp)
            form.line_bonds_connected_root_2.setText(comp)
        if paths is not None: 
            form.textedit_bonds_active_filters_2.append('Shortest paths from '+paths[0]+' to '+paths[1])
            form.line_bonds_path_root_2.setText(paths[0])
            form.line_bonds_path_goal_2.setText(paths[1])
        if avg_least is not None: 
            form.textedit_bonds_active_filters_2.append('Minimal bonds paths from '+avg_least[0]+' to '+avg_least[1])
            form.line_bonds_path_root_2.setText(avg_least[0])
            form.line_bonds_path_goal_2.setText(avg_least[1])
        if segm is not None: 
            if segm[1] is not None:
                form.textedit_bonds_active_filters_2.append('Between segments '+segm[0]+' and '+str(segm[1]))
                form.line_bonds_filter_sega_2.setText(segm[0])
                form.line_bonds_filter_segb_2.setText(segm[1])
            else:
                form.textedit_bonds_active_filters_2.append('Between segment '+segm[0]+' and all')
                form.line_bonds_filter_sega_2.setText(segm[0])
        if resn is not None: 
            if resn[1] is not None:
                form.textedit_bonds_active_filters_2.append('Between resnames '+resn[0]+' and '+str(resn[1]))
                form.line_bonds_filter_resa_2.setText(resn[0])
                form.line_bonds_filter_resb_2.setText(resn[1])
            else:
                form.textedit_bonds_active_filters_2.append('Between resname '+resn[0]+' and all')
                form.line_bonds_filter_resa_2.setText(resn[0])
        
        
    form_plotting.button_wire_draw_water_in_wire.clicked.connect(draw_wire_water_in_wire)
    #form_plotting.button_wires_remove_last.clicked.connect(remove_wire_joint_occupancy)
    form_plotting.button_wires_add_wwa.clicked.connect(compare_wire_joint_occupancy)
    form_plotting.button_wires_draw_occupancy.clicked.connect(draw_wire_occupancy)
    form_plotting.button_wires_compare_occupancy.clicked.connect(compare_wire_occupancy)
    form_plotting.button_wires_draw_bonds_per_residue.clicked.connect(draw_wire_per_residue)
    form_plotting.button_wires_compare_bonds_per_residue.clicked.connect(compare_wire_per_residue)
    form_plotting.button_wires_draw_joint_barcode.clicked.connect(draw_wire_joint_occupancy)
    form_plotting.button_wires_draw_bonds_timeseries.clicked.connect(draw_wire_sum_timeseries)
    form_plotting.button_wires_compare_bonds_timeseries.clicked.connect(compare_wire_sum_timeseries)
    form_plotting.button_wires_draw_bonds_per_time_interval.clicked.connect(draw_wire_per_time_interval)
    form_plotting.button_wires_compare_bonds_per_time_interval.clicked.connect(compare_wire_per_time_interval)
    form_plotting.button_wires_draw_heatmap.clicked.connect(draw_wire_heatmap)
    form_plotting.checkbox_wires_heatmap_residuewise.stateChanged.connect(toggle_wire_heatmap_lines)
    form_plotting.checkbox_wires_occupancy_segnames.stateChanged.connect(toggle_wire_occupancy_segments)
    form_plotting.checkbox_wires_timeseries_segnames.stateChanged.connect(toggle_wire_timeseries_segments)
    form.button_wire_open_plotting.clicked.connect(load_plotting_wires)
    form.button_wire_draw_graph.clicked.connect(draw_wire_graph)
    #form.button_wire_projection.clicked.connect(compute_wire_projection)
    #form.button_wire_average.clicked.connect(compute_average)
    #form.button_wire_joint.clicked.connect(compute_wires_joint)
    form.button_bonds_filter_2.clicked.connect(filter_wires)
    form.button_wires_clear.clicked.connect(clear_all_filters_wires)
    form.button_wire_structure.clicked.connect(browse_wire_structure)
    form.button_wire_trajectories.clicked.connect(browse_wire_trajectories)
    form.button_wire_init.clicked.connect(init_wire_analysis)
    form.button_wire_set.clicked.connect(set_wires)
    form.button_wire_save.clicked.connect(save_wires_to_file)
    form.button_wire_restore.clicked.connect(restore_wires_from_file)
    #form.button_wire_graph_save.clicked.connect(save_wire_graph)
    form.checkBox_wires_angle.stateChanged.connect(toggle_wire_angle)
    form.pushButton_wires_advanced.clicked.connect(load_advanced_wires)
    form_plotting.pushButton.clicked.connect(save_plot)
    form_plotting.button_plotting_save_data.clicked.connect(save_data)
    form.pushButton_wires_save_data.clicked.connect(save_wires_data)
    form_graph_wires.pushButton_bonds_save_data.clicked.connect(save_wires_data_advanced)
    form_graph_wires.pushButton_bonds_compare_graphs.clicked.connect(draw_wire_graph_compare)
    form_graph_wires.pushButton_draw_graph.clicked.connect(draw_wire_graph_advanced)
    write_license()
    
    dialog.show()