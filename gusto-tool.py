#!/usr/bin/env python


import os
from PySide import QtGui, QtCore


class FileBrowser(QtGui.QWidget):

    class Signals(QtCore.QObject):
        files_loaded = QtCore.Signal(list, list)

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
        tree_view.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)

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



class DatasetBrowser(QtGui.QWidget):

    class Signals(QtCore.QObject):
        dset_selected = QtCore.Signal(str)

    def __init__(self):
        super(DatasetBrowser, self).__init__()

        remove_button = QtGui.QPushButton('Remove')
        remove_button.clicked.connect(self.handle_remove_button)

        model = QtGui.QStandardItemModel()
        model.itemChanged.connect(self.handle_item_changed)

        list_view = QtGui.QListView()
        list_view.setModel(model)

        selection_model = list_view.selectionModel()
        selection_model.currentChanged.connect(self.handle_current_changed)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(list_view)
        layout.addWidget(remove_button)

        self.setLayout(layout)
        self.model = model
        self.signals = self.Signals()


    def get_checked_indexes(self):
        checked_indexes = [ ]
        i = 0
        while self.model.item(i):
            if self.model.item(i).checkState():
                checked_indexes.append(i)
            i += 1
        checked_indexes.sort(reverse=True)
        return checked_indexes


    def get_items(self):
        items = [ ]
        i = 0
        while self.model.item(i):
            items.append(self.model.item(i))
            i += 1
        return items


    def handle_remove_button(self):
        for i in self.get_checked_indexes():
            self.model.removeRows(i, 1)


    def handle_item_changed(self, item):
        print "items were changed", item


    def handle_current_changed(self, current, previous):
        item = self.model.item(current.row())
        if item:
            self.signals.dset_selected.emit(item.filename)


    def load_files(self, files=[], dirs=[]):
        filenames = [item.filename for item in self.get_items()]

        for f in files:
            if f in filenames: continue
            item = QtGui.QStandardItem(os.path.basename(f))
            item.filename = f
            item.setToolTip(f)
            item.setCheckable(True)
            item.setSelectable(True)
            self.model.appendRow(item)

        for directory in dirs:
            files = [f for f in os.listdir(directory) if f.endswith('.h5')]
            self.load_files(files)



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



class PlottingArea(QtGui.QWidget):

    def __init__(self):
        super(PlottingArea, self).__init__()

        from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
        from matplotlib.figure import Figure
        from matplotlib.axes import Axes
        from matplotlib.lines import Line2D

        reset_button = QtGui.QPushButton('Reset')
        reset_button.clicked.connect(self.reset_plot)

        fig = Figure()
        ax1 = Axes(fig, [0.1, 0.1, 0.8, 0.8])
        fig.add_axes(ax1)

        canvas = FigureCanvas(fig)
        navbar = NavigationToolbar(canvas, self)
        radio_group = RadioButtonGroup(['X', 'R', 'z', 'u0', 'u1', 'u2', 'u3', 'dg', 'pg', '-alfv', '-fast'])
        radio_group.signals.selection_changed.connect(self.change_variable)

        layout = QtGui.QGridLayout()
        layout.addWidget(canvas, 0, 1)
        layout.addWidget(navbar, 1, 1)
        layout.addWidget(reset_button, 3, 1)
        layout.addWidget(radio_group, 0, 0, QtCore.Qt.AlignTop)
        self.setLayout(layout)

        self.fig = fig
        self.axes = ax1
        self.canvas = canvas
        self.navbar = navbar

        self.variable = 'u0'
        self.dset_filename = None
        self.reset_plot()


    def plot_dset(self, dset_filename):
        self.dset_filename = dset_filename
        self.replot()


    def change_variable(self, new_variable):
        self.variable = new_variable
        self.replot()


    def reset_plot(self):
        from matplotlib.lines import Line2D
        self.axes.cla()
        self.line = Line2D([], [])
        self.axes.add_line(self.line)
        self.replot()


    def replot(self):
        import gusto_dataset

        if not self.dset_filename: return
        if not self.variable: return

        dset = gusto_dataset.GustoDataset(self.dset_filename)
        x = dset.get_cell_variable('R')
        y = dset.get_cell_variable(self.variable)

        self.line.set_data(x, y)
        self.axes.relim(visible_only=True)
        self.axes.autoscale_view()
        self.axes.set_title(self.dset_filename)
        self.axes.set_xlabel('R')
        self.axes.set_ylabel(self.variable)
        self.canvas.draw()

        dset.close()



class MainWindow(QtGui.QWidget):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.resize(1200, 800)
        self.setWindowTitle("Gusto tool")

        #file_browser = FileBrowser()
        dset_browser = DatasetBrowser()
        plotting_area = PlottingArea()

        #file_browser.signals.files_loaded.connect(dset_browser.load_files)
        dset_browser.signals.dset_selected.connect(plotting_area.plot_dset)
        dset_browser.load_files(dirs=[os.getcwd()])

        layout = QtGui.QHBoxLayout()
        #layout.addWidget(file_browser)
        layout.addWidget(dset_browser)
        layout.addWidget(plotting_area)
        self.setLayout(layout)



if __name__ == '__main__':
    import sys
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    window.raise_()
    app.exec_()
