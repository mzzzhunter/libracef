import sys
import pandas as pd
import matplotlib.pyplot as plt
from PyQt5.QtGui import QIcon, QKeySequence
from PyQt5.QtWidgets import QApplication, QMainWindow, QTableView, QVBoxLayout, QWidget, QFileDialog, QMenu, QMessageBox, QHBoxLayout, QAction, QSplitter, QInputDialog, QHeaderView, QDialog, QLabel, QDoubleSpinBox, QDialogButtonBox, QLineEdit, QPushButton, QListWidget, QComboBox, QTableWidget, QTableWidgetItem
from PyQt5.QtCore import Qt, QAbstractTableModel
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from nist_util import search_nist_for_spectrum
from export_util import create_jcamp_library, create_mslibrary_xml
from cef_util import combine_cef_results
import ast
import os
from pathlib import Path

class PandasModel(QAbstractTableModel):
    def __init__(self, data):
        QAbstractTableModel.__init__(self)
        self._data = data
        self._decimal_places = {}

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole or role == Qt.EditRole:
                value = self._data.iloc[index.row(), index.column()]
                if isinstance(value, float):
                    decimal_places = self._decimal_places.get(index.column(), 2)
                    return f"{value:.{decimal_places}f}"
                return str(value)
        return None

    def setData(self, index, value, role):
        if role == Qt.EditRole:
            try:
                self._data.iloc[index.row(), index.column()] = value
                self.dataChanged.emit(index, index)
                return True
            except ValueError:
                return False
        return False

    def headerData(self, col, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self._data.columns[col]
        return None

    def flags(self, index):
        return Qt.ItemIsSelectable | Qt.ItemIsEnabled | Qt.ItemIsEditable

    def set_decimal_places(self, column, places):
        self._decimal_places[column] = places
        self.layoutChanged.emit()

class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

class CEFImportDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("CEF Import Settings")
        layout = QVBoxLayout(self)

        # Add directory selection
        dir_layout = QHBoxLayout()
        self.dir_input = QLineEdit()
        dir_layout.addWidget(QLabel("CEF Directory:"))
        dir_layout.addWidget(self.dir_input)
        browse_button = QPushButton("Browse")
        browse_button.clicked.connect(self.browse_directory)
        dir_layout.addWidget(browse_button)
        layout.addLayout(dir_layout)

        self.rt_tolerance_input = QDoubleSpinBox()
        self.rt_tolerance_input.setRange(0, 1)
        self.rt_tolerance_input.setSingleStep(0.01)
        self.rt_tolerance_input.setValue(0.1)
        layout.addWidget(QLabel("RT Tolerance:"))
        layout.addWidget(self.rt_tolerance_input)

        self.group_similarity_threshold_input = QDoubleSpinBox()
        self.group_similarity_threshold_input.setRange(0, 1)
        self.group_similarity_threshold_input.setSingleStep(0.01)
        self.group_similarity_threshold_input.setValue(0.9)
        layout.addWidget(QLabel("Group Similarity Threshold:"))
        layout.addWidget(self.group_similarity_threshold_input)

        # Add Exclude Element input
        self.exclude_element_input = QLineEdit()
        self.exclude_element_input.setText("Si")
        layout.addWidget(QLabel("Exclude Element (separate by ','):"))
        layout.addWidget(self.exclude_element_input)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def browse_directory(self):
        directory = QFileDialog.getExistingDirectory(self, "Select CEF Directory")
        if directory:
            self.dir_input.setText(directory)

    def get_values(self):
        exclude_elements = [elem.strip() for elem in self.exclude_element_input.text().split(',')]
        return (self.dir_input.text(), self.rt_tolerance_input.value(),
                self.group_similarity_threshold_input.value(), exclude_elements)
# Create a dialog box for selecting column and order
class SortDialog(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.setWindowTitle("Sort By")

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.column_combo = QComboBox()
        self.column_combo.addItems(self.parent().df.columns)
        self.layout.addWidget(self.column_combo)

        self.order_combo = QComboBox()
        self.order_combo.addItems(['Ascending', 'Descending'])
        self.layout.addWidget(self.order_combo)

        self.button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        self.layout.addWidget(self.button_box)

    def get_values(self):
        return self.column_combo.currentText(), self.order_combo.currentText()

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.df = None
        self.bar_width_percent = 0.5  # Default bar width as 0.5% of canvas width
        self.setWindowIcon(QIcon('libracef_icon.jpg'))
        self.initUI()
        self.curent_csv_file = None
        self.nist_path = self.find_nist_ms_search_default_paths()
        self.undo_stack = []

    def initUI(self):
        self.setWindowTitle('LibraCEF - Build MS library from CEFs')
        self.setGeometry(100, 100, 1200, 600)

        main_layout = QHBoxLayout()
        
        left_layout = QVBoxLayout()
        
        self.table = QTableView()
        self.table.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self.show_context_menu)
        self.table.horizontalHeader().setContextMenuPolicy(Qt.CustomContextMenu)
        self.table.horizontalHeader().customContextMenuRequested.connect(self.show_header_context_menu)
        self.table.setSelectionMode(QTableView.ExtendedSelection)  # Enable multiple row selection
        self.table.clicked.connect(self.plot_current_spectrum)
        left_layout.addWidget(self.table)

        left_widget = QWidget()
        left_widget.setLayout(left_layout)
        
        self.canvas = MplCanvas(self, width=5, height=4, dpi=100)
        
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(left_widget)
        splitter.addWidget(self.canvas)
        splitter.setSizes([2*self.width()//3, self.width()//3])
        
        main_layout.addWidget(splitter)

        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

        # Add File menu
        menubar = self.menuBar()
        file_menu = menubar.addMenu('File')
        import_csv_action = QAction('Import CSV', self)
        import_csv_action.triggered.connect(self.import_csv)
        file_menu.addAction(import_csv_action)
        
        import_cef_action = QAction('Import CEF results', self)
        import_cef_action.triggered.connect(self.import_cef)
        file_menu.addAction(import_cef_action)

        # Add a separator
        file_menu.addSeparator()

        self.import_save_action = QAction("Save", self)
        self.import_save_action.triggered.connect(self.save_csv)
        file_menu.addAction(self.import_save_action)
        self.import_save_action.setEnabled(False)
                
        self.export_csv_action = QAction('Save as', self)
        self.export_csv_action.triggered.connect(self.export_to_csv)
        file_menu.addAction(self.export_csv_action)
        self.export_csv_action.setEnabled(False)

        # Add Export menu
        export_menu = menubar.addMenu('Export')
        
        export_jcamp_action = QAction('Export to JCAMP Library', self)
        export_jcamp_action.triggered.connect(self.export_to_jcamp)
        export_menu.addAction(export_jcamp_action)

        export_mslibrary_action = QAction('Export to MSLibrary XML', self)
        export_mslibrary_action.triggered.connect(self.export_to_mslibrary)
        export_menu.addAction(export_mslibrary_action)

        # Add Settings menu
        settings_menu = menubar.addMenu('Settings')
        set_bar_width_action = QAction('Set Bar Width', self)
        set_bar_width_action.triggered.connect(self.set_bar_width)
        settings_menu.addAction(set_bar_width_action)

        # Add action to set NIST path
        set_nist_path_action = QAction('Set NIST MS Search Path', self)
        set_nist_path_action.triggered.connect(self.set_nist_path)
        settings_menu.addAction(set_nist_path_action)

        # Add Edit menu
        edit_menu = menubar.addMenu('Edit')
        ## Add an undo action
        undo_action = QAction("Undo", self)
        undo_action.triggered.connect(self.undo_changes)
        edit_menu.addAction(undo_action)

        ## Add a shortcut to the undo action
        undo_action.setShortcut(QKeySequence.Undo)

        table_menu = edit_menu.addMenu('Table')
        add_column_action = QAction('Add Column', self)
        add_column_action.triggered.connect(self.add_column)
        table_menu.addAction(add_column_action)
        
        remove_column_action = QAction('Remove Column', self)
        remove_column_action.triggered.connect(self.remove_column)
        table_menu.addAction(remove_column_action)
        
        rearrange_column_action = QAction('Rearrange Columns', self)
        rearrange_column_action.triggered.connect(self.rearrange_columns)
        table_menu.addAction(rearrange_column_action)

        # Add the sort by action
        sort_by_action = QAction('Sort By', self)
        sort_by_action.triggered.connect(self.sort_by_column)
        table_menu.addAction(sort_by_action)

        # Add Help menu
        help_menu = menubar.addMenu('Help')
        about_action = QAction('About', self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)
        
        
        #add status bar
        self.statusBar = self.statusBar()

        # Add a permanent label to the status bar
        self.status_label = QLabel("DataFrame is empty")
        self.statusBar.addWidget(self.status_label)

    def save_csv(self):
        if self.df is not None and self.curent_csv_file is not None:
            try:
                self.df.to_csv(self.curent_csv_file, index=False)
                # QMessageBox.information(self, "Save Successful", f"Data saved to {self.current_file}")
            except Exception as e:
                QMessageBox.warning(self, "Save Error", f"Failed to save file: {str(e)}")
                print(f"Error saving file: {str(e)}")
        elif self.df is not None:
            QMessageBox.warning(self, "Save Error", "No file loaded. Please load a file first.")
        else:
            QMessageBox.warning(self, "Save Error", "No data to save.")

    def import_csv(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Open CSV File", "", "CSV Files (*.csv)")

        if file_name:
            self.df = pd.read_csv(file_name, keep_default_na=False)
            self.update_table()
            # Update the window title with the file name
            self.curent_csv_file = file_name
            self.set_window_title(file_name)
            self.import_save_action.setEnabled(True)
            self.export_csv_action.setEnabled(True)

    def import_cef(self):
        dialog = CEFImportDialog(self)
        if dialog.exec_():
            directory, rt_tolerance, group_similarity_threshold, exclude_elements = dialog.get_values()
            if directory:
                try:
                    self.df = combine_cef_results(directory, rt_tolerance=rt_tolerance, 
                                                group_similarity_threshold=group_similarity_threshold,
                                                exclude_elements=exclude_elements)
                    self.update_table()
                    # self.import_save_action.setEnabled(True)
                    self.export_csv_action.setEnabled(True)
                except:
                    QMessageBox.warning(self, "Import Error", "Error when import CEF files.")

            else:
                QMessageBox.warning(self, "Import Error", "Please select a valid directory.")

    def update_table(self):
        if self.df is not None:
            model = PandasModel(self.df)
            self.table.setModel(model)
            # show the number of rows and columns in the status bar
            rows, cols = self.df.shape
            self.status_label.setText(f"DataFrame updated: {rows} rows, {cols} columns")

    def show_context_menu(self, pos):
        context = QMenu(self)
        # show_spectrum_action = context.addAction("Show Spectrum")
        search_nist_action = context.addAction("Search NIST")
        insert_row_action = context.addAction("Insert Row")
        delete_rows_action = context.addAction("Delete Selected Rows")
        modify_ms_spectrum_action = context.addAction("Modify MS spectrum")       
        action = context.exec_(self.table.mapToGlobal(pos))
        
        index = self.table.indexAt(pos)
        if action == delete_rows_action:
            self.delete_selected_rows()
        elif action == search_nist_action:            
            self.search_nist(index.row())
        elif action == insert_row_action:
            self.insert_row(index.row())
        elif action == modify_ms_spectrum_action:
            self.modify_ms_spectrum(index.row())

    def show_header_context_menu(self, pos):
        context = QMenu(self)
        set_decimal_places_action = context.addAction("Set Decimal Places")
        action = context.exec_(self.table.horizontalHeader().mapToGlobal(pos))

        if action == set_decimal_places_action:
            column = self.table.horizontalHeader().logicalIndexAt(pos)
            self.set_decimal_places(column)

    def set_decimal_places(self, column):
        places, ok = QInputDialog.getInt(self, "Set Decimal Places", 
                                         "Enter number of decimal places:",
                                         value=2, min=0, max=10)
        if ok:
            model = self.table.model()
            if isinstance(model, PandasModel):
                model.set_decimal_places(column, places)

    def plot_current_spectrum(self):
        if self.df is None:
            return
        index = self.table.currentIndex()
        self.plot_spectrum(index.row())

    def plot_spectrum(self, row):
        if self.df is None:
            return
        label_offset = 10
        ms_peaks = self.df.iloc[row]['MS_Peaks']
        
        # Check if MS_Peaks is a string and convert it to a list if necessary
        if isinstance(ms_peaks, str):
            try:
                ms_peaks = ast.literal_eval(ms_peaks)
            except:
                QMessageBox.warning(self, "Plot Error", "Unable to parse MS_Peaks data.")
                return
        
        if not isinstance(ms_peaks, list):
            QMessageBox.warning(self, "Plot Error", "MS_Peaks data is not in the expected format.")
            return
        
        mz, intensity = zip(*ms_peaks)
        
        self.canvas.axes.clear()
        
        # Calculate bar width as set percentage of the canvas width
        canvas_width = self.canvas.figure.get_figwidth() * self.canvas.figure.get_dpi()
        bar_width = self.bar_width_percent / 100 * (max(mz) - min(mz))
        # print(mz, intensity)
        bars = self.canvas.axes.bar(mz, intensity, width=bar_width)
        self.canvas.axes.set_xlabel('m/z')
        self.canvas.axes.set_ylabel('Intensity')
        self.canvas.axes.set_title(f"{self.df.iloc[row]['Chemical_Name']}")
        self.canvas.axes.set_ylim(0, 1200)
        
        # Add x value (m/z) on top of each bar
        for mz, intensity in zip(mz, intensity):
            self.canvas.axes.text(mz + bar_width/2., intensity + label_offset, f'{round(mz)}', ha='center', va='bottom', rotation=90, fontsize=8)      
        self.canvas.draw()

    def insert_row(self, selected_row):
        # Append Undo stack
        undo_state = self.df.copy()
        self.undo_stack.append(undo_state)
        
        new_row =  self.df.iloc[selected_row].copy()
        new_row['Chemical_Name'] = "New Compound"
        self.df = pd.concat([self.df.iloc[:selected_row], pd.DataFrame([new_row]), self.df.iloc[selected_row:]]).reset_index(drop=True)
        self.update_table()


    def delete_selected_rows(self):
        if self.df is None:
            return

        selected_rows = sorted(set(index.row() for index in self.table.selectedIndexes()), reverse=True)
        
        if not selected_rows:
            QMessageBox.warning(self, "Delete Error", "No rows selected.")
            return

        reply = QMessageBox.question(self, 'Confirm Deletion',
                                     f"Are you sure you want to delete {len(selected_rows)} selected row(s)?",
                                     QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            undo_state = self.df.copy()
            self.undo_stack.append(undo_state)
            for row in selected_rows:
                self.df = self.df.drop(self.df.index[row])
            self.df = self.df.reset_index(drop=True)
            self.update_table()

    def undo_changes(self):
        if self.undo_stack:
            self.df = self.undo_stack.pop()
            self.update_table()

    def search_nist(self, row):
        if self.df is None:
            return
        if self.nist_path:
            current_path = Path(os.path.abspath(os.getcwd()))
            spec_data_path = current_path / 'spectrum.txt' # add spectrum file to the path
            df_row = self.df.iloc[row]
            
            search_nist_for_spectrum(df_row, spec_data_path=str(spec_data_path), nist_path=str(self.nist_path))
        else:            
            QMessageBox.warning(self, 'NIST MS Search path not defined', 
                             'Please define the NIST MS Search path in the Settings before searching.')
            

    def export_to_csv(self):
        if self.df is None:
            QMessageBox.warning(self, "Export Error", "No data to export.")
            return

        file_name, _ = QFileDialog.getSaveFileName(self, "Save CSV File", "", "CSV Files (*.csv)")
        if file_name:
            self.df.to_csv(file_name, index=False)
            self.curent_csv_file = file_name
            self.set_window_title(file_name)
            self.import_save_action.setEnabled(True)
            QMessageBox.information(self, "Export Successful", f"Data exported to {file_name}")

    def export_to_jcamp(self):
        if self.df is None:
            QMessageBox.warning(self, "Export Error", "No data to export.")
            return

        file_name, _ = QFileDialog.getSaveFileName(self, "Save JCAMP Library", "", "JCAMP Files (*.jdx)")
        if file_name:
            jcamp_library = create_jcamp_library(self.df)
            with open(file_name, 'w') as f:
                for entry in jcamp_library:
                    f.write(entry + '\n\n')
            QMessageBox.information(self, "Export Successful", f"JCAMP library exported to {file_name}")
    def export_to_mslibrary(self):
        if self.df is None:
            QMessageBox.warning(self, "Export Error", "No data to export.")
            return

        file_name, _ = QFileDialog.getSaveFileName(self, "Save MSLibrary XML", "", "XML Files (*.mslibrary.xml)")
        if file_name:
            if not file_name.endswith('.mslibrary.xml'):
                file_name = os.path.splitext(file_name)[0] + '.mslibrary.xml'
            mslibrary_xml = create_mslibrary_xml(self.df, file_name)
            with open(file_name, 'w') as f:
                f.write(mslibrary_xml)

    def set_bar_width(self):
        width, ok = QInputDialog.getDouble(self, "Set Bar Width", 
                                           "Enter bar width as percentage of canvas width (0.1-10):",
                                           self.bar_width_percent, 0.1, 10, 2)
        if ok:
            self.bar_width_percent = width
            # If there's a current plot, update it
            current_row = self.table.currentIndex().row()
            if current_row >= 0:
                self.plot_spectrum(current_row)

    def set_nist_path(self):
    # Open a file dialog to select the NIST folder
        nist_path = QFileDialog.getExistingDirectory(self, 'Select NIST MS Search Folder')
        if nist_path:
            # Save the NIST path to a config file or a variable
            # Use pathlib to make the path compatible for different OS
            nist_path = Path(nist_path)

            # Check if nistms$.exe exist in the folder
            file_path = nist_path / 'nistms$.exe'
            if file_path.is_file():
                self.nist_path = nist_path
            else:
                QMessageBox.warning(self, "NIST MS Search Path Error", "Cannot find nistms$.exe")

    def add_column(self):
        if not self._check_df_exists():
            return
        column_name, ok = QInputDialog.getText(self, "Add Column", "Enter new column name:")
        if ok and column_name:
            self.df[column_name] = None
            self.update_table()

    def remove_column(self):
        if not self._check_df_exists():
            return
        column_names = list(self.df.columns)
        column_name, ok = QInputDialog.getItem(self, "Remove Column", "Select a column to remove:", column_names, 0, False)
        if ok and column_name:
            self.df = self.df.drop(columns=[column_name])
            self.update_table()

    def rearrange_columns(self):
        if self.df is not None:
            column_names = list(self.df.columns)
            
            rearrange_dialog = QDialog(self)
            rearrange_dialog.setWindowTitle("Rearrange Columns")
            rearrange_dialog.setGeometry(100, 100, 400, 400)

            vertical_layout = QVBoxLayout()
            rearrange_dialog.setLayout(vertical_layout)

            column_list = QListWidget()
            for name in column_names:
                column_list.addItem(name)

            column_list_item = column_list.currentItem()

            def move_up():
                nonlocal column_list_item
                current_row = column_list.currentRow()
                if current_row > 0:
                    current_text = column_list.currentItem().text()
                    column_list.takeItem(current_row)
                    column_list.insertItem(current_row - 1, current_text)
                    column_list.setCurrentRow(current_row - 1)
                    column_list_item = column_list.currentItem()

            def move_down():
                nonlocal column_list_item
                current_row = column_list.currentRow()
                if current_row < column_list.count() - 1:
                    current_text = column_list.currentItem().text()
                    column_list.takeItem(current_row)
                    column_list.insertItem(current_row + 1, current_text)
                    column_list.setCurrentRow(current_row + 1)
                    column_list_item = column_list.currentItem()

            def save_changes():
                nonlocal column_names
                new_column_order = []
                for row in range(column_list.count()):
                    new_column_order.append(column_list.item(row).text())
                self.df = self.df[new_column_order]
                self.update_table()
                rearrange_dialog.close()

            horizontal_layout1 = QHBoxLayout()
            move_up_button = QPushButton("Move Up")
            move_up_button.clicked.connect(move_up)
            horizontal_layout1.addWidget(move_up_button)

            move_down_button = QPushButton("Move Down")
            move_down_button.clicked.connect(move_down)
            horizontal_layout1.addWidget(move_down_button)

            save_button = QPushButton("Save Changes")
            save_button.clicked.connect(save_changes)
            vertical_layout.addWidget(column_list)
            vertical_layout.addLayout(horizontal_layout1)
            vertical_layout.addWidget(save_button)

            rearrange_dialog.exec_()
        else:
            QMessageBox.warning(self, "Rearrange Columns Error", "Please import or create a data frame first.")

    # Create a function to sort the table by column and order
    def sort_by_column(self):
        if self.df is not None:
            dialog = SortDialog(self)
            if dialog.exec_():
                column_name, order = dialog.get_values()
                
                self.df= self.df.sort_values(by=column_name, ascending=(order == 'Ascending'))
                self.update_table()
                # # Example: show the number of rows and columns in the status bar
                # rows, cols = self.df.shape
                # self.status_label.setText(f"DataFrame updated: {rows} rows, {cols} columns")


    def modify_ms_spectrum(self, selected_row):
        
        if selected_row < 0:
            QMessageBox.warning(self, "No Selection", "Please select a row to modify.")
            return

        ms_peaks = self.df.at[selected_row, 'MS_Peaks']

        # Check if MS_Peaks is a string and convert it to a list if necessary 
        if isinstance(ms_peaks, str):
            try:
                ms_peaks = ast.literal_eval(ms_peaks)
            except:
                QMessageBox.warning(self, "Plot Error", "Unable to parse MS_Peaks data.")
                return     
               
        if not isinstance(ms_peaks, list):
            ms_peaks = []

        dialog = QDialog(self)
        dialog.setWindowTitle("Modify MS Spectrum")
        layout = QVBoxLayout(dialog)

        table = QTableWidget(dialog)
        table.setColumnCount(2)
        table.setHorizontalHeaderLabels(["m/z", "Intensity"])
        table.setRowCount(len(ms_peaks))

        for i, peak in enumerate(ms_peaks):
            table.setItem(i, 0, QTableWidgetItem(str(peak[0])))
            table.setItem(i, 1, QTableWidgetItem(str(peak[1])))

        add_button = QPushButton("Add Row", dialog)
        remove_button = QPushButton("Remove Row", dialog)
        save_button = QPushButton("Save", dialog)

        def add_row():
            table.insertRow(table.rowCount())

        def remove_row():
            if table.rowCount() > 0:
                table.removeRow(table.currentRow())

        def save_changes():
            new_ms_peaks = []
            for row in range(table.rowCount()):
                mz_item = table.item(row, 0)
                intensity_item = table.item(row, 1)
                if mz_item and intensity_item:
                    mz = float(mz_item.text())
                    intensity = float(intensity_item.text())
                    new_ms_peaks.append((mz, intensity))
            self.df.at[selected_row, 'MS_Peaks'] = new_ms_peaks
            self.plot_spectrum(selected_row)
            dialog.accept()

        add_button.clicked.connect(add_row)
        remove_button.clicked.connect(remove_row)
        save_button.clicked.connect(save_changes)

        layout.addWidget(table)
        layout.addWidget(add_button)
        layout.addWidget(remove_button)
        layout.addWidget(save_button)

        dialog.setLayout(layout)
        dialog.exec()

    def show_about(self):
        QMessageBox.about(self, "About", "LibraCEF v0.1\n\nAuthor: Tony Chen\n\nhttps://github.com/mzzzhunter/libracef\n\nLicense: MIT License")

    def set_window_title(self, file_name):
        base_title = "LibraCEF"
        file_name = os.path.basename(file_name)
        self.setWindowTitle(f"{base_title} - {file_name}")

    def find_nist_ms_search_default_paths(self):

        windir = os.getenv("windir")
        if not windir:
            return None

        win_ini_path = os.path.join(windir, "win.ini")
        if not os.path.exists(win_ini_path):
            return None
        
        path32_value = None
        with open(win_ini_path, 'r') as file:
            for line in file:
                if line.strip().lower().startswith("path32="): 
                    path32_value = line.split("=", 1)[1].strip()
                    
        return path32_value
    def _check_df_exists(self):
        if self.df is None:
            QMessageBox.warning(self, "Error", "Please import or create a data frame first.")
            return False
        return True

def main():
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()























