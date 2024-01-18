#!/usr/bin/env python
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import QFileDialog
import sys
import sqlite3
import os

palm_dir = os.getcwd()
palm_bin = palm_dir + '/trunk/SCRIPTS'


class MyTable(QtWidgets.QDialog):

    class MyTableWidgetItem(QtWidgets.QTableWidgetItem):
        def __init__(self, text, sortKey):
            QtWidgets.QTableWidgetItem.__init__(self, text, QtWidgets.QTableWidgetItem.UserType)
            self.sortKey = sortKey

        def __lt__(self, other):
            return self.sortKey < other.sortKey

    def __init__(self):
        super(MyTable, self).__init__()
        uic.loadUi(palm_bin + '/palm_gf_files/palm_gf_table.ui', self)
        self.center_window()
        self.load_limits()
        self.exec_button.clicked.connect(self.filter_data)
        # connects labels of table to sorting subroutine
        self.nx_table.clicked.connect(lambda: self.sort_table(str("nx")))
        self.ny_table.clicked.connect(lambda: self.sort_table(str("ny")))
        self.nz_table.clicked.connect(lambda: self.sort_table(str("nz")))
        self.npex_table.clicked.connect(lambda: self.sort_table(str("npex")))
        self.npey_table.clicked.connect(lambda: self.sort_table(str("npey")))
        self.npexnpey_table.clicked.connect(lambda: self.sort_table(str("pxpy")))
        self.np_table.clicked.connect(lambda: self.sort_table(str("np")))
        self.ngpts_table.clicked.connect(lambda: self.sort_table(str("ngpts")))
        self.nxpex_table.clicked.connect(lambda: self.sort_table(str("nxnpex")))
        self.nypey_table.clicked.connect(lambda: self.sort_table(str("nynpey")))
        # connects main buttons to subroutines
        self.load_all_button.clicked.connect(self.load_all)
        self.reset_button.clicked.connect(self.load_limits)
        self.save_button.clicked.connect(self.save)

    def center_window(self):
        # centers window in the middle of screen
        frame = self.frameGeometry()
        center = QtWidgets.QDesktopWidget().availableGeometry().center()
        frame.moveCenter(center)
        self.move(frame.topLeft())

    def load_limits(self):
        # reads out the .db and sets limits of results as Min/Max
        con = sqlite3.connect('palm_gf_data.db')
        c = con.cursor()
        c.execute("SELECT * FROM " + 'grid_limits')
        Min = c.fetchone()
        Max = c.fetchone()

        self.nx_min.setValue(Min[0])
        self.nx_min.setMinimum(Min[0])
        self.nx_max.setValue(Max[0])
        self.nx_max.setMaximum(Max[0])

        self.ny_min.setValue(Min[1])
        self.ny_min.setMinimum(Min[1])
        self.ny_max.setValue(Max[1])
        self.ny_max.setMaximum(Max[1])

        self.nz_min.setValue(Min[2])
        self.nz_min.setMinimum(Min[2])
        self.nz_max.setValue(Max[2])
        self.nz_max.setMaximum(Max[2])

        self.npex_min.setValue(Min[3])
        self.npex_min.setMinimum(Min[3])
        self.npex_max.setValue(Max[3])
        self.npex_max.setMaximum(Max[3])

        self.npey_min.setValue(Min[4])
        self.npey_min.setMinimum(Min[4])
        self.npey_max.setValue(Max[4])
        self.npey_max.setMaximum(Max[4])

        self.npxnpy_min.setValue(Min[5])
        self.npxnpy_min.setMinimum(Min[5])
        self.npxnpy_max.setValue(Max[5])
        self.npxnpy_max.setMaximum(Max[5])

        self.np_min.setValue(Min[6])
        self.np_min.setMinimum(Min[6])
        self.np_max.setValue(Max[6])
        self.np_max.setMaximum(Max[6])
        # try statement, as ngpts could be too huge to handle for PyQt5 (>10^8)
        try:
            self.ngpts_min.setValue(Min[7])
            self.ngpts_min.setMinimum(Min[7])
            self.ngpts_max.setValue(Max[7])
            self.ngpts_max.setMaximum(Max[7])

        except OverflowError:
            self.ngpts_min.setEnabled(False)
            self.ngpts_max.setEnabled(False)
            self.ngpts_table.setEnabled(False)
            self.ngpts_label.setStyleSheet('color: grey')
            self.exec_button.setEnabled(False)
            self.status_label.setText('ngpts overflow')

        self.nxpex_min.setValue(Min[8])
        self.nxpex_min.setMinimum(Min[8])
        self.nxpex_max.setValue(Max[8])
        self.nxpex_max.setMaximum(Max[8])

        self.nypey_min.setValue(Min[9])
        self.nypey_min.setMinimum(Min[9])
        self.nypey_max.setValue(Max[9])
        self.nypey_max.setMaximum(Max[9])

        con.commit()
        c.close()
        con.close()

    def load_all(self):
        # shows all results found in .db file
        con = sqlite3.connect('palm_gf_data.db')
        c = con.cursor()
        c.execute("SELECT * FROM " + 'grid_current')
        data = c.fetchall()

        self.tableWidget.setRowCount(len(data))
        row_cnt = -1
        for i in range(0, len(data)):
            line = data[i]
            row_cnt += 1
            for j in range(0, 10):
                value = line[j]

                if j == 7:
                    self.tableWidget.setItem(row_cnt, j, self.MyTableWidgetItem(
                        str("%.1e" % value), i))

                else:
                    self.tableWidget.setItem(row_cnt, j, self.MyTableWidgetItem(
                        str(value), i))

        c.close()
        con.close()
        self.tableWidget.setRowCount(row_cnt + 1)
        self.status_label.setText(str(row_cnt + 1) + ' results loaded')

    def filter_data(self):
        # filters data according to set Min/Max. Reads out .db to find all fitting data
        con = sqlite3.connect('palm_gf_data.db')
        c = con.cursor()
        c.execute("SELECT * FROM " + 'grid_current')
        data = c.fetchall()

        self.tableWidget.setRowCount(len(data))

        row_cnt = -1
        for i in range(0, len(data)):
            line = data[i]
            if self.nx_min.value() <= line[0] <= self.nx_max.value():
                if self.ny_min.value() <= line[1] <= self.ny_max.value():
                    if self.nz_min.value() <= line[2] <= self.nz_max.value():
                        if self.npex_min.value() <= line[3] <= self.npex_max.value():
                            if self.npey_min.value() <= line[4] <= self.npey_max.value():
                                if self.npxnpy_min.value() <= line[5] <= self.npxnpy_max.value():
                                    if self.np_min.value() <= line[6] <= self.np_max.value():
                                        if self.ngpts_min.value() <= line[7] <= self.ngpts_max.value():
                                            if self.nxpex_min.value() <= line[8] <= self.nxpex_max.value():
                                                if self.nypey_min.value() <= line[9] <= self.nypey_max.value():
                                                    row_cnt += 1
                                                    for j in range(0, 10):
                                                        value = line[j]

                                                        if j == 7:
                                                            self.tableWidget.setItem(row_cnt, j, self.MyTableWidgetItem(
                                                                str("%.1e" % value), i))
                                                        else:
                                                            self.tableWidget.setItem(row_cnt, j, self.MyTableWidgetItem(
                                                                str(value), i))

        c.close()
        con.close()
        self.tableWidget.setRowCount(row_cnt + 1)
        self.status_label.setText(str(row_cnt + 1) + ' results found')

    def sort_table(self, column):
        # with argument as chosen column, read not from table but from .db and reloads sorted data into table
        fnx_mn = self.nx_min.value()
        fnx_mx = self.nx_max.value()
        fny_mn = self.ny_min.value()
        fny_mx = self.ny_max.value()
        fnz_mn = self.nz_min.value()
        fnz_mx = self.nz_max.value()
        fnpex_mn = self.npex_min.value()
        fnpex_mx = self.npex_max.value()
        fnpey_mn = self.npex_min.value()
        fnpey_mx = self.npey_max.value()
        fnpxnpy_mn = self.npxnpy_min.value()
        fnpxnpy_mx = self.npxnpy_max.value()
        fnp_mn = self.np_min.value()
        fnp_mx = self.np_max.value()
        fngpts_mn = self.ngpts_min.value()
        fngpts_mx = self.ngpts_max.value()
        nxpex_mn = self.nxpex_min.value()
        nxpex_mx = self.nxpex_max.value()
        nypey_mn = self.nypey_min.value()
        nypey_mx = self.nypey_max.value()

        if column == 'nx':
            sorted_col = 'nx'
            if self.nx_table.isChecked() is True:
                order = " DESC"
            else:
                order = " ASC"
            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("ny"):
            sorted_col = "ny"
            if self.ny_table.isChecked() is True:
                order = " DESC"
            else:
                order = " ASC"

            self.nx_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("nz"):
            sorted_col = "nz"
            if self.nz_table.isChecked() is True:
                order = " DESC"
            else:
                order = " ASC"

            self.nx_table.setChecked(False)
            self.ny_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("npex"):
            sorted_col = "npex"
            if self.npex_table.isChecked() is True:
                order = " DESC"
            else:
                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("npey"):
            sorted_col = "npey"
            if self.npey_table.isChecked() is True:
                order = " DESC"
            else:
                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("pxpy"):
            sorted_col = "pxpy"
            if self.npexnpey_table.isChecked() is True:
                order = " DESC"
            else:
                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.np_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("np"):
            sorted_col = "np"
            if self.np_table.isChecked() is True:
                order = " DESC"
            else:
                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("ngpts"):
            sorted_col = "ngpts"
            if self.ngpts_table.isChecked() is True:
                order = " DESC"
            else:
                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("nxnpex"):
            sorted_col = "nxnpex"
            if self.nxpex_table.isChecked() is True:
                order = " DESC"
            else:
                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("nynpey"):
            sorted_col = "nynpey"
            if self.nypey_table.isChecked() is True:
                order = " DESC"
            else:
                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)

        conn = sqlite3.connect('palm_gf_data.db')
        c = conn.cursor()
        c.execute("SELECT * FROM grid_current  WHERE nx <= " + str(fnx_mx) + " AND nx >= " + str(
            fnx_mn) + " AND ny <= " + str(fny_mx) + " AND ny >= " + str(fny_mn) + " AND nz <= " + str(fnz_mx) +
                  " AND nz >= " + str(fnz_mn) + " AND npex <= " + str(fnpex_mx) + " AND npex >= " +
                  str(fnpex_mn) + " AND npey <= " + str(fnpey_mx) + " AND npey >= " + str(fnpey_mn) + " AND "
                                                                                                      "pxpy <= " + str(
            fnpxnpy_mx) + " AND pxpy >= " + str(fnpxnpy_mn) + " AND np <= " + str(fnp_mx) + " AND np >= " + str(
            fnp_mn) + " AND ngpts <= " + str(fngpts_mx) + " AND ngpts >= " + str(fngpts_mn) +
                  " AND nxnpex <= " + str(nxpex_mx) + " AND nxnpex >= " + str(nxpex_mn) + " AND nynpey <= " + str(
            nypey_mx) + " AND nynpey >= " + str(nypey_mn) +
                  " ORDER BY " + str(sorted_col) + str(order))

        sorted = c.fetchall()
        c.close()
        conn.close()
        self.tableWidget.setRowCount(len(sorted))
        for row_indx in range(0, len(sorted)):
            for col_indx in range(0, 10):
                row = sorted[row_indx]
                value = row[col_indx]
                if col_indx == 7:
                    self.tableWidget.setItem(row_indx, col_indx, self.MyTableWidgetItem(str("%.1e" % value), col_indx))
                else:
                    self.tableWidget.setItem(row_indx, col_indx, self.MyTableWidgetItem(str(value), col_indx))

        self.status_label.setText(str(len(sorted)) + ' results found')

    def save_dialog(self):
        # calling native PyQt5 Dialog to get desired save path
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self, "Save .txt File", ".txt",
                                                  "All Files (*.txt);;Text Files (*.txt)", options=options)
        if fileName:
            return fileName

    def save(self):
        # saves the file
        # try because QFileDialog could be closed without choosing path
        try:
            if len(self.tableWidget.selectedItems()) != 0:
                location = self.save_dialog()
                data2safe = self.tableWidget.selectedItems()
                if location is not None:
                    with open(location, 'w') as file:
                        file.write('nx  ' + 'ny  ' + 'nz  ' + 'npex  ' + 'npey' + "\n")
                        for row in range(0, int(len(data2safe)), 10):
                            file.write(str(data2safe[0+row].text()) + '  ' + str(data2safe[1+row].text()) + '  ' +
                                       str(data2safe[2+row].text()) + '  ' + str(data2safe[3+row].text()) + '  ' +
                                       str(data2safe[4+row].text()) + "\n")
                    file.close()
            else:
                self.status_label.setText('Select to save')
        except AttributeError:
            pass


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MyTable()
    window.setWindowTitle('palm_gf Table')
    window.show()
    sys.exit(app.exec_())
