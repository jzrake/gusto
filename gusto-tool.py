#!/usr/bin/env python


import os
from PySide import QtGui, QtCore


class FileBrowser(QtGui.QWidget):

    class Signals(QtCore.QObject):
        files_loaded = QtCore.Signal(list, list)
        file_changed = QtCore.Signal(str)

    def __init__(self, parent=None):
        super(FileBrowser, self).__init__(parent)

        model = QtGui.QFileSystemModel()
        model.setReadOnly(True)
        model.setNameFilters(["*.h5"])
        model.setNameFilterDisables(False)
        model.directoryLoaded.connect(self.resize_tree_view)
        model.setRootPath(os.getcwd())

        tree_view = QtGui.QTreeView()
        tree_view.setModel(model)
        tree_view.setRootIndex(model.index(os.getcwd()))
        tree_view.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        tree_view.hideColumn(1)
        tree_view.hideColumn(2)
        tree_view.hideColumn(3)

        selection_model = tree_view.selectionModel()
        selection_model.currentChanged.connect(self.handle_current_changed)

        load_button = QtGui.QPushButton('Load')
        load_button.clicked.connect(self.handle_load_button)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(tree_view)
        layout.addWidget(load_button)

        self.setLayout(layout)
        self.model = model
        self.tree_view = tree_view
        self.signals = self.Signals()


    def handle_load_button(self):
        dirs = [ ]
        files = [ ]
        for index in self.tree_view.selectedIndexes():
            if index.column() != 0: continue
            filename = self.model.filePath(index)
            if os.path.isdir(filename):
                dirs.append(filename)
            else:
                files.append(filename)
        self.signals.files_loaded.emit(files, dirs)


    def resize_tree_view(self):
        for i in range(4):
            self.tree_view.resizeColumnToContents(i)


    def handle_current_changed(self, index):
        filename = self.model.filePath(index)
        self.signals.file_changed.emit(filename)



class RadioButtonGroup(QtGui.QWidget):

    class Signals(QtCore.QObject):
        selection_changed = QtCore.Signal(str)

    def __init__(self, labels, parent=None):
        super(RadioButtonGroup, self).__init__(parent)

        layout = QtGui.QVBoxLayout()
        buttons = [ ]

        for label in labels:
            button = QtGui.QRadioButton(label)
            button.toggled.connect(self.handle_toggled)
            layout.addWidget(button)
            buttons.append(button)

        self.setLayout(layout)
        self.buttons = buttons
        self.signals = self.Signals()

    def handle_toggled(self, enabled):
        if not enabled: return
        for button in self.buttons:
            if button.isChecked():
                self.signals.selection_changed.emit(button.text())



class CheckboxGroup(QtGui.QWidget):

    class Signals(QtCore.QObject):
        state_changed = QtCore.Signal(list)

    def __init__(self, labels, parent=None):
        super(CheckboxGroup, self).__init__(parent)

        layout = QtGui.QVBoxLayout()
        buttons = [ ]

        for label in labels:
            button = QtGui.QCheckBox(label)
            button.stateChanged.connect(self.handle_state_changed)
            layout.addWidget(button)
            buttons.append(button)

        self.setLayout(layout)
        self.buttons = buttons
        self.signals = self.Signals()

    def handle_state_changed(self, state):
        self.signals.state_changed.emit([b.text() for b in self.buttons
                                         if b.checkState()])



class PlottingArea(QtGui.QWidget):

    def __init__(self):
        super(PlottingArea, self).__init__()

        from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
        from matplotlib.figure import Figure
        from matplotlib.axes import Axes

        reset_button = QtGui.QPushButton('Reset')
        reset_button.clicked.connect(self.reset_lines)

        fig = Figure()
        ax1 = Axes(fig, [0.1, 0.1, 0.8, 0.8])
        fig.add_axes(ax1)

        canvas = FigureCanvas(fig)
        navbar = NavigationToolbar(canvas, self)

        all_variables = ['X', 'R', 'z',
                         'u0', 'u2', 'u3', 'b2', 'b3',
                         'dg', 'pg', 'entropy', 'sigma',
                         '-slow', '-alfv', '-fast',
                         'flux:mass',
                         'flux:energy',
                         'flux:angmom']

        var_group = CheckboxGroup(all_variables)
        var_group.signals.state_changed.connect(self.set_active_variables)

        opt_group = CheckboxGroup(['log x', 'log y'])
        opt_group.signals.state_changed.connect(self.set_bool_options)

        all_markers = ['', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd']
        mrk_combo = QtGui.QComboBox()

        for mrk in all_markers:
            mrk_combo.addItem(mrk)
        mrk_combo.currentIndexChanged.connect(lambda i: self.set_marker(mrk_combo.itemText(i)))

        lw_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        lw_slider.setValue(20)
        lw_slider.setMinimum(5)
        lw_slider.valueChanged.connect(self.set_line_width)

        layout = QtGui.QGridLayout()
        layout.addWidget(canvas, 0, 1, 4, 1)
        layout.addWidget(navbar, 4, 1)
        layout.addWidget(var_group, 0, 0, QtCore.Qt.AlignTop)
        layout.addWidget(opt_group, 1, 0, QtCore.Qt.AlignBottom)
        layout.addWidget(mrk_combo, 2, 0, QtCore.Qt.AlignBottom)
        layout.addWidget(lw_slider, 3, 0, QtCore.Qt.AlignBottom)
        layout.addWidget(reset_button, 4, 0)
        self.setLayout(layout)

        self.fig = fig
        self.axes = ax1
        self.canvas = canvas
        self.navbar = navbar

        self.lines = { }
        self.active_variables = [ ]
        self.all_variables = all_variables
        self.active_bools = [ ]
        self.dset_filename = None
        self.reset_lines()


    def set_dset_filename(self, dset_filename):
        if dset_filename.endswith('.h5'):
            self.dset_filename = dset_filename
        else:
            self.dset_filename = None
        self.reload_data()
        self.draw()


    def set_active_variables(self, variables):
        self.active_variables = variables
        self.reload_data()
        self.draw()


    def set_bool_options(self, options):
        self.active_bools = options
        self.draw()


    def set_line_width(self, lw):
        for line in self.lines.itervalues():
            line.set_linewidth(lw / 10.)
        self.canvas.draw()


    def set_marker(self, marker):
        for line in self.lines.itervalues():
            line.set_marker(marker)
        self.canvas.draw()


    def reset_lines(self):
        from matplotlib.lines import Line2D
        self.axes.cla()
        for key in self.all_variables:
            self.lines[key] = Line2D([], [])
            self.axes.add_line(self.lines[key])
        self.reload_data()
        self.draw()


    def reload_data(self):
        import gusto_dataset

        if self.dset_filename == None:
            for line in self.lines.itervalues():
                line.set_visible(False)
            return

        dset = gusto_dataset.GustoDataset(self.dset_filename)
        for key in self.lines:
            line = self.lines[key]
            if key in self.active_variables:
                x = dset.get_cell_variable('R')
                y = dset.get_cell_variable(key)
                line.set_data(x, y)
                line.set_visible(True)
                line.set_color(self.make_color_from_string_hash(key))
            else:
                line.set_visible(False)
        dset.close()


    def draw(self):
        ylabel = ''
        title = os.path.basename(self.dset_filename) if self.dset_filename else ''

        self.axes.set_xscale('log' if 'log x' in self.active_bools else 'linear')
        self.axes.set_yscale('log' if 'log y' in self.active_bools else 'linear')
        self.axes.relim(visible_only=True)
        self.axes.autoscale_view()
        self.axes.set_title(title)
        self.axes.set_xlabel('R')
        self.axes.set_ylabel(ylabel)
        self.canvas.draw()




    def make_color_from_string_hash(self, key):
        r = (hash(key) >> 0*8) % 256
        g = (hash(key) >> 1*8) % 256
        b = (hash(key) >> 2*8) % 256
        return [r/256., g/256., b/256.]


class MainWindow(QtGui.QWidget):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.resize(1200, 800)
        self.setWindowTitle("Gusto tool")

        file_browser = FileBrowser()
        plotting_area = PlottingArea()

        file_browser.signals.file_changed.connect(plotting_area.set_dset_filename)

        layout = QtGui.QHBoxLayout()
        layout.addWidget(file_browser)
        layout.addWidget(plotting_area)
        self.setLayout(layout)



if __name__ == '__main__':
    import sys
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    window.raise_()
    app.exec_()
