#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------#
# This file is part of the PALM model system.
#
# PALM is free software: you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 1997-2021  Leibniz Universitaet Hannover
#--------------------------------------------------------------------------------#
#
# Description:
# ------------
# Graphical user interface for the palmrun script.
# @author Felix Gaschler
# @author Björn Maronga (maronga@muk.uni-hannover.de)
#
# Instructions:
# -------------
#
#------------------------------------------------------------------------------!
import sys
import subprocess
from PyQt5 import QtCore, QtGui, uic, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5.QtCore import QProcess, pyqtSlot, QTimer
from time import strftime
import os
import shutil

global version
version = "PALM 23.10-rc.1"
global update_frequency
update_frequency = 5*60000

devnull = open(os.devnull, 'w')
install_bin_dir = os.path.dirname(__file__)
palmrungui_sdir = os.path.dirname(os.path.realpath(__file__))
palm_dir = os.path.dirname(install_bin_dir)
palm_source = os.path.join(palmrungui_sdir, '../../palm/model/src')
job_dir = os.path.join(palm_dir, 'JOBS')

if not os.path.exists(job_dir):
   os.mkdir(job_dir)

# returns the human readable file size of filename
def file_size(filename, suffix='B'):
   size_raw = os.stat(filename).st_size
   for unit in [' ',' K',' M',' G',' T']:
      if abs(size_raw) < 1024.0:
         return "%i%s%s" % (size_raw, unit, suffix)
      size_raw /= 1024.0
   return "%.1f%s%s" % (num, 'Y', suffix)

def ilen(iterable):
    return reduce(lambda sum, element: sum + 1, iterable, 0)

palmrunline = ""
set_list = []

Ui_MainWindow = uic.loadUiType("%s/palmrungui_files/mainwindow.ui" % (palmrungui_sdir))[0]
Ui_helpDialog = uic.loadUiType("%s/palmrungui_files/help.ui" % (palmrungui_sdir))[0]
Ui_aboutDialog = uic.loadUiType("%s/palmrungui_files/about.ui" % (palmrungui_sdir))[0]

class HelpDialog(QtWidgets.QDialog,Ui_helpDialog):
    def __init__(self, parent=None):
        super(HelpDialog,self).__init__()
        self.setupUi(self)
        
class AboutDialog(QtWidgets.QDialog,Ui_aboutDialog):
    def __init__(self, parent=None):
        super(AboutDialog,self).__init__()
        self.setupUi(self)        

class Mainwindow(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super(Mainwindow, self).__init__()
        self.setupUi(self)

        self.palm_logo.setPixmap(QtGui.QPixmap(palmrungui_sdir + "/palmrungui_files/logo.png"))
        self.groupBox.findChild(QtWidgets.QLabel,"label_version").setText(version)
        self.recent_jobs(50)
        self.load_jobs()



        # look up configuration files and add to combo box
        
        self.group_execution.findChild(QtWidgets.QComboBox,"combo_configuration").addItem("")
        
        for files in os.listdir(palm_dir):
           if files.startswith(".palm.config"):
              tmpstring = list(filter(None,files.split(".palm.config.")))[0]
              self.group_execution.findChild(QtWidgets.QComboBox,"combo_configuration").addItem(tmpstring)

        commandline = self.groupBox.findChild(QtWidgets.QLineEdit,"commandline")
        commandline.setText("")
        
        self.tabWidget.setCurrentIndex(0)  
        

        filename = "%s/.palmrungui.default" % (palm_dir)  
        if os.path.exists(filename):
           file = open(filename, "r")
           if ( file is not None ):
              # File opened successfully
              palmrunline = file.readline()
              #if neue zeile zeichen
              palmrunline = palmrunline[:len(palmrunline)-1]
              file.close()

           # In case a palmrunline was found, load it to mainwindow
           if ( palmrunline != ""):
              palmrunline = palmrunline[17:]
              commandline.setText(palmrunline)
              self.setup_gui(palmrunline)

        QtWidgets.QApplication.processEvents()

        # Start refresh timer. On timeout perform update

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_all)
        self.timer.setSingleShot(False)
        self.timer.start(update_frequency)

        # The timetimer counts the time since last update
        self.timetimer= QtCore.QElapsedTimer()
        self.timetimer.start()

        # The labeltimer induces the update of the remaining time in the UI

        self.labeltimer = QtCore.QTimer(self)
        self.labeltimer.timeout.connect(self.UpdatePush)
        self.labeltimer.setSingleShot(False)

        # Update in the UI will be performed at each 1/10 of the update interval
        self.labeltimer.start(int(update_frequency/10))
        self.push_update.setText("Update (" + str(int(update_frequency/1000/60)) + " min)")

    # starts xterm with palmrun commandline
    #######################################
    def startpalmrun(self):
        palmrunline = str(self.groupBox.findChild(QtWidgets.QLineEdit,"commandline").text())
        palmrunline_save = palmrunline

        # Disable the main window
        self.tabWidget.setEnabled(False)
        self.groupBox.findChild(QtWidgets.QPushButton,"button_start").setEnabled(False)
        self.groupBox.findChild(QtWidgets.QLineEdit,"line_tag").setEnabled(False)
        self.groupBox.findChild(QtWidgets.QPushButton,"button_start").setText("wait...")
        self.groupBox.findChild(QtWidgets.QLineEdit,"commandline").setText("Executing palmrun in xterm...")

        # Wait until all commands have been executed (ugly) ?
        #for i in range(0,21):
        #    qApp->processEvents()       

        # Start xterm as QProcess
        palmrun = QProcess()
        palmrun.setProcessChannelMode(QProcess.MergedChannels) # mergedChannels
        palmrun.setWorkingDirectory(palm_dir)

        geomet = self.frameGeometry()

        posx = geomet.x()+geomet.width()
        posy = geomet.y()

        if ( os.path.isfile("palmrungui.log") ):
           os.remove("palmrungui.log")

        s = " -title \"Executing palmrun...\" -fa \"Monospace\" -l -lf \"palmrungui.log\" -fs 11 -geometry \"80x38+%d+%d\" -e \"" % (posx,posy)
        palmrunline = "%s%s ;echo -n '--> Press enter to continue...';read yesno\"</dev/stdin" % (s, palmrunline.replace("\"","\'"))

        mString = "xterm %s" % (palmrunline)
        palmrun.start(mString)

        if( palmrun.waitForStarted() is False ):
            return

        # Wait until palmrun has finished or wait for 200 minutes
        palmrun.waitForFinished(3600000)        

        # Jobs has been submitted or aborted. Continuing...
        # Save the palmrun command to history file
        filename = "%s/.palmrungui.history" % (palm_dir)
        tmstamp = strftime("%Y/%m/%d %H:%M")
        tag = str(self.groupBox.findChild(QtWidgets.QLineEdit,"line_tag").text())
    
        file = open(filename,"a")
        s = "%s %s (%s)\n" % (tmstamp,palmrunline_save,tag)
        file.write(s)
        file.close()              
        
        # Enable main window again
        self.tabWidget.setEnabled(True)
        self.groupBox.findChild(QtWidgets.QPushButton,"button_start").setEnabled(True)
        self.groupBox.findChild(QtWidgets.QPushButton,"button_start").setText("palmrun")
        self.groupBox.findChild(QtWidgets.QLineEdit,"commandline").setText(palmrunline_save)
        self.groupBox.findChild(QtWidgets.QLineEdit,"line_tag").setEnabled(True)
        
        # Reload recent jobs
        self.recent_jobs(50)
        
        
       
    # Update the label
    def UpdatePush(self):
        remaining_time = int((update_frequency - self.timetimer.elapsed()) / 1000 / 60)
        self.push_update.setText("Update (" + str(remaining_time) + " min)")


    # starts xterm with palmbuild commandline
    #######################################
    def startpalmbuild(self):
        palmbuildline = 'palmbuild -c ' + str(self.group_execution.findChild(QtWidgets.QComboBox, "combo_configuration").currentText())

        full_build = self.group_execution.findChild(QtWidgets.QCheckBox,"check_rebuild").checkState()

        if ( full_build == 2 ):
           for filename in os.listdir(palm_source):
              os.utime(palm_source + "/" + filename, None)

        palmrunline_save = str(self.groupBox.findChild(QtWidgets.QLineEdit,"commandline").text())

        # Disable the main window
        palmbuildline_save = palmbuildline
        self.tabWidget.setEnabled(False)
        self.groupBox.findChild(QtWidgets.QPushButton,"button_start").setEnabled(False)
        self.groupBox.findChild(QtWidgets.QLineEdit,"line_tag").setEnabled(False)
        self.groupBox.findChild(QtWidgets.QPushButton,"button_start").setText("wait...")
        self.group_execution.findChild(QtWidgets.QPushButton,"button_palmbuild").setText("wait...")
        self.groupBox.findChild(QtWidgets.QLineEdit,"commandline").setText("Executing palmbuild in xterm...")
        
        # Wait until all commands have been executed (ugly) ?
        #for i in range(0,21):
        #    qApp->processEvents()       
        
        # Start xterm as QProcess
        palmbuild = QProcess()
        palmbuild.setProcessChannelMode(QProcess.MergedChannels) # mergedChannels
        palmbuild.setWorkingDirectory(palm_dir)
    
        geomet = self.frameGeometry()
        
        posx = geomet.x()+geomet.width()
        posy = geomet.y()
      
        s = " -title \"Executing palmbuild...\" -fa \"Monospace\" -fs 11 -geometry \"80x38+%d+%d\" -e \"" % (posx,posy)
        palmbuildline = "%s%s;echo -n '--> Press enter to continue...';read yesno\"</dev/stdin" % (s, palmbuildline.replace("\"","\'"))
        
        mString = "xterm %s" % (palmbuildline)
        palmbuild.start(mString)
    
        if( palmbuild.waitForStarted() is False ):
            return
        
        # Wait until palmbuild has finished or wait for 200 minutes
        palmbuild.waitForFinished(3600000)        
        
          
        
        # Enable main window again
        self.tabWidget.setEnabled(True)
        self.groupBox.findChild(QtWidgets.QPushButton,"button_start").setEnabled(True)
        self.groupBox.findChild(QtWidgets.QLineEdit,"line_tag").setEnabled(True)
        self.groupBox.findChild(QtWidgets.QPushButton,"button_start").setText("palmrun")
        self.groupBox.findChild(QtWidgets.QLineEdit,"commandline").setText(palmrunline_save)
        self.group_execution.findChild(QtWidgets.QPushButton,"button_palmbuild").setText("(re-)build")
        # Reload recent jobs
        self.recent_jobs(50)

    # loads recent jobs 
    ######################################
    def recent_jobs(self, number):
        fileDir = "%s/.palmrungui.history" % (palm_dir)
        if os.path.exists(fileDir):
            pass
        else:
            return

        file = open(fileDir,"r")
        history = file.readlines()
        tmphistory = list()
        tmptag = list()
        file.close()
        j = 0

        list_jobname = self.group_history.findChild(QtWidgets.QListWidget,"list_jobname")
        list_jobname.clear()

        # Read history entries and append to recent job list
        i=len(history)-1
        count = 0
        while i>=0 and count<number:
            timestamp = history[i][:16]
            listitem = history[i][17:len(history[i])-1]
            tagitem = listitem
            configitem = listitem
            matchitems = list_jobname.findItems(listitem, QtCore.Qt.MatchExactly)

            if ( len(matchitems) == 0 ):
                listitem = list(filter(None,listitem.split(" -r")))[1]
                listitem = listitem.strip()
                listitem = list(filter(None,listitem.split(" ")))[0]
                listitem = listitem.replace("\"","") 
                list_jobname.addItem(listitem)
                
                configitem = list(filter(None,configitem.split(" -c")))
                if ( len(configitem) < 2):
                   count = count +1
                   continue
                else:
                   configitem = configitem[1].strip()
                configitem = list(filter(None,configitem.split(" ")))[0]
                configitem = configitem.replace("\"","") 

                s = "%s: %s (%s)" % (timestamp,listitem,configitem)
                tmphistory.append(s)
                
                
                tagitem = list(filter(None,tagitem.split(" (")))[1]
                tagitem = tagitem.replace(")","")  
                if ( len(tagitem) == 0 ): 
                   s = "Tag: empty"
                else:
                   s = "Tag: %s" % (tagitem)                
                tmptag.append(s)
                count = count +1

                j = j+1
                
            if ( j == number ):
                break
            i = i-1
            
        # Send to list
        list_jobname.clear()
        
        i=0
        while i<len(tmphistory):
            list_jobname.addItem(tmphistory[i])
            list_jobname.item(i).setToolTip(tmptag[i])
            i = i+1
            

    # Enables coupled settings
    ###############################
    def enable_coupled(self):
        coupledState = self.group_execution.findChild(QtWidgets.QComboBox, "drop_job").currentText()
        group = self.group_execution.findChild(QtWidgets.QGroupBox, "group_coupled")
        
        if (coupledState == "Restart run (coupled atmosphere ocean)" or coupledState == "Initial run (coupled atmosphere ocean)"):
            self.group_coupled.setEnabled(True)
        else: 
            self.group_coupled.setEnabled(False)
            

    # select a job via click from list
    #################################
    def choosejob_list(self):
        #  Get selected item from list
        list_jobname = self.group_history.findChild(QtWidgets.QListWidget,"list_jobname")
        filename = str(list_jobname.currentItem().text())

        timestamp = filename[:16]
        jobname = filename[18:].split(" ")[0]

        itemint = list_jobname.currentRow()

        # Reload list
        self.recent_jobs(50)
        self.load_jobs()
    
        # Set selected item to jobname
        list_jobname.item(itemint).setSelected(True)
    
        # Set selected item to jobname in "available jobs" list
        item2int = self.list_jobs.findItems(jobname,QtCore.Qt.MatchCaseSensitive)
        
        if ( item2int != [] ):
           self.list_jobs.setCurrentItem(item2int[0])
        self.update_input()
           
    
        fileDir = "%s/.palmrungui.history" % (palm_dir)
        file = open(fileDir,"r")
        history = file.readlines()
        tmphistory = list()
        file.close()        
    
        i = len(history)
        while i>=1:
            listitem = history[i-1][17:len(history[i-1])-1]
            listitem = list(filter(None,listitem.split(" -r")))[1]
            listitem = listitem.strip()
            listitem = list(filter(None,listitem.split(" ")))[0]
            listitem = listitem.replace("\"","") 

            # Select command line with correct timestamp and jobname
            if (history[i-1][:16] == timestamp and listitem == jobname):
                palmrunline = history[i-1]
                palmrunline = palmrunline[17:]
                palmrunline = palmrunline[:len(palmrunline)-1]
                palmrunline = list(filter(None,palmrunline.split("(")))[0]
                self.groupBox.findChild(QtWidgets.QLineEdit,"commandline").setText(palmrunline)
                
              
                tag = history[i-1].split('\n')[0]
                tag = list(filter(None,tag.split("(")))[1]
                tag = tag.replace(")","")

                self.groupBox.findChild(QtWidgets.QLineEdit,"line_tag").setText(tag)
                self.setup_gui(palmrunline)
                self.list_jobname.item(itemint).setSelected(True)    
                
                return
            i = i-1
                
    # Change run identifer (initial, restart, coupled...)
    ######################################################
    def change_rc_list(self):
        drop_job = self.group_execution.findChild(QtWidgets.QComboBox,"drop_job").currentText()
        self.change_commandline("a","")
    
        # Enable PE distribution for atmosphere/ocean
        if ( drop_job == "Restart run (coupled atmosphere ocean)" or drop_job == "Initial run (coupled atmosphere ocean)"):
            drop_atmos = self.group_coupled.findChild(QtWidgets.QLineEdit,"line_PE_atmos").text()
            drop_ocean = self.group_coupled.findChild(QtWidgets.QLineEdit,"line_PE_ocean").text()
            s = "%s %s" % (drop_atmos,drop_ocean)
            self.change_commandline("Y",s)
    
    
        if ( drop_job == "Restart run" or drop_job == "Restart run (coupled atmosphere ocean)"):
           self.change_commandline("r","")

        # Check of ocean runs
        else:
            self.delete_commandline("Y")
            if (drop_job == "Precursor run (ocean)"):
                self.activate_flag("y")
            else:
                self.deactivate_flag("y")
            
    # changes commandline depending on parameters
    ##########################################################      
    def change_commandline(self, id_str, fwt_str):
        fwt_str = str(fwt_str) 
        initialize = False
        palmrunline = str(self.groupBox.findChild(QtWidgets.QLineEdit,"commandline").text())
        s = " -%s " % (id_str)
        splitline = list(filter(None,palmrunline.split(s)))

        if ( len(splitline) == 0 ):
                splitline.append("palmrun")
                splitline.append(" ")
                initialize = True

        elif ( len(splitline) == 1 ):
            splitline.append(" ")
        
        param = splitline[1].split("-")
    
        # Change in parameter "r" (jobname)
        if (id_str == "r"):
            filename = str(self.group_execution.findChild(QtWidgets.QLineEdit,"line_jobname").text())
            s = filename.split("JOBS/") 
            param[0] = s[len(s)-1]
     
            if ( initialize == True ):#and self.group_runcontrol.isEnabled() == True ):
                self.group_execution.setEnabled(True)
                self.group_execution.findChild(QtWidgets.QComboBox,"drop_job").setEnabled(True)
                self.group_advanced.setEnabled(True)
                self.groupBox.findChild(QtWidgets.QPushButton,"button_start").setEnabled(True)
                self.menuBar.findChild(QtWidgets.QMenu,"menuStart").actions()[3].setEnabled(True)
                self.groupBox.findChild(QtWidgets.QLineEdit,"line_tag").setEnabled(True)
            
            elif ( param[0] == ""):            
                self.group_execution.setEnabled(False)
                self.group_execution.findChild(QtWidgets.QComboBox,"drop_job").setEnabled(False)
                self.groupBox.findChild(QtWidgets.QPushButton,"button_start").setEnabled(False)
                self.menuBar.findChild(QtWidgets.QMenu,"menuStart").actions()[3].setEnabled(False)
                self.group_advanced.setEnabled(False)
                self.groupBox.findChild(QtWidgets.QLineEdit,"line_tag").setEnabled(True)
                self.delete_commandline("r")
                self.change_commandline("a","remove")
                return 1 
            
            else:

                self.group_execution.setEnabled(True)
                self.drop_job.setEnabled(True)
                self.group_advanced.setEnabled(True)
                   
        # Change in parameter "a" (run control list)
        elif (id_str == "a"): 
    
            drop_job = self.group_execution.findChild(QtWidgets.QComboBox,"drop_job").currentText()
            rc_flag = "#"
    
            if (drop_job == "Initial run"):
                rc_flag = "#"
            
            elif (drop_job == "Restart run"):            
                rc_flag = "r"
            
            elif (drop_job == "Precursor run (atmosphere)"):            
                rc_flag = "#"
            
            elif (drop_job == "Precursor run (ocean)"):            
                rc_flag = "o#"

            elif (drop_job == "Initial run (coupled atmosphere ocean)"):            
                rc_flag = "#"
            
            elif (drop_job == "Restart run (coupled atmosphere ocean)"):            
                rc_flag = "r"
            
            param[0] = "\"d3%s" % (rc_flag)
    
    
            if (drop_job == "Restart run (coupled atmosphere ocean)" or drop_job == "Initial run (coupled atmosphere ocean)"):
                if (rc_flag == "#"):
                   rc_flag = "o#"
                else:
                   rc_flag = "or"

                param[0] = "%s d3%s" % (param[0],rc_flag)
    
            status_restarts = self.group_execution.findChild(QtWidgets.QCheckBox,"check_restarts").checkState()

            if (status_restarts == 2):
                param[0]="%s restart" % (param[0])


            status_cycfill = self.group_execution.findChild(QtWidgets.QCheckBox,"check_cycfill").checkState()
    
            if (status_cycfill == 2):            
                param[0]="%s rec" % (param[0])

            status_svf = self.group_execution.findChild(QtWidgets.QCheckBox,"check_svf").checkState()

            if (status_svf == 2):            
               param[0]="%s svfin" % (param[0])
 
            status_svf_2 = self.group_execution.findChild(QtWidgets.QCheckBox,"check_svf_2").checkState()

            if (status_svf_2 == 2):
               param[0]="%s svfout" % (param[0])

            param[0]="%s\"" % (param[0])

            if ( fwt_str == "remove"):            
                self.delete_commandline(id_str)
                return 1
            
            else:            
                self.button_start.setEnabled(True)
                self.action_save.setEnabled(True)
                self.groupBox.findChild(QtWidgets.QLineEdit,"line_tag").setEnabled(True)

        # Change in any other parameter
        else:
            if ( fwt_str != ""):            
                param[0] = "\"%s\"" % (fwt_str)
            
            else:            
                self.delete_commandline(id_str)
                return 1
                        
        # Join the new palmrunline
        splitline[1]= " -".join(param)
                
                
        s = " -%s " % (id_str)
        newpalmrunline = s.join(splitline)
    
        # Print the new palmrunline to mainwindow
        newpalmrunline = newpalmrunline.replace("  "," ")
        self.groupBox.findChild(QtWidgets.QLineEdit,"commandline").setText(newpalmrunline)

    # change lineinput depending on sender
    ###################################################################################
    def change_lineinput(self):
        if ( self.sender() == self.group_execution.findChild(QtWidgets.QComboBox, "combo_configuration") ):
            tmptext = self.group_execution.findChild(QtWidgets.QComboBox, "combo_configuration").currentText()
            if not tmptext:
               self.change_commandline("c"," ")
               self.group_execution.findChild(QtWidgets.QPushButton,'button_palmbuild').setEnabled(False)
            else:
               self.change_commandline("c",tmptext)
               self.group_execution.findChild(QtWidgets.QPushButton,'button_palmbuild').setEnabled(True)

        elif ( self.sender() == self.group_execution.findChild(QtWidgets.QLineEdit,"line_q")):
            tmptext = self.group_execution.findChild(QtWidgets.QLineEdit,"line_q").text()
            self.change_commandline("q",tmptext)
        
        elif ( self.sender() == self.group_execution.findChild(QtWidgets.QLineEdit,"line_account")):
            tmptext = self.group_execution.findChild(QtWidgets.QLineEdit,"line_account").text()
            self.change_commandline("A",tmptext)

        elif ( self.sender() == self.group_execution.findChild(QtWidgets.QLineEdit,"line_wait")):
            tmptext = self.group_execution.findChild(QtWidgets.QLineEdit,"line_wait").text()
            self.change_commandline("W",tmptext)

        elif ( self.sender() ==  self.group_execution.findChild(QtWidgets.QLineEdit,"line_pe")):
            tmptext = self.group_execution.findChild(QtWidgets.QLineEdit,"line_pe").text()
            self.change_commandline("X",tmptext)
        
        elif ( self.sender() == self.group_execution.findChild(QtWidgets.QLineEdit,"line_tpn")):
            tmptext = self.group_execution.findChild(QtWidgets.QLineEdit,"line_tpn").text()
            self.change_commandline("T",tmptext)
                
        elif ( self.sender() == self.group_execution.findChild(QtWidgets.QLineEdit,"line_time")):
            tmptext = self.group_execution.findChild(QtWidgets.QLineEdit,"line_time").text()
            self.change_commandline("t",tmptext)
        
        elif ( self.sender() == self.group_advanced.findChild(QtWidgets.QLineEdit,"line_M")):
            tmptext = self.group_advanced.findChild(QtWidgets.QLineEdit,"line_M").text()
            self.change_commandline("M",tmptext)
        
        elif ( self.sender() == self.group_advanced.findChild(QtWidgets.QLineEdit,"line_m")):
            tmptext = self.group_advanced.findChild(QtWidgets.QLineEdit,"line_m").text()
            self.change_commandline("m",tmptext)
        
        elif ( self.sender() == self.group_advanced.findChild(QtWidgets.QLineEdit,"line_D")):
            tmptext = self.group_advanced.findChild(QtWidgets.QLineEdit,"line_D").text()
            self.change_commandline("D",tmptext)
        
        elif ( self.sender() == self.group_advanced.findChild(QtWidgets.QLineEdit,"line_c")):
            tmptext = self.group_advanced.findChild(QtWidgets.QLineEdit,"line_c").text()
            if ( tmptext == ".palmrungui.config"):
                tmptext = ""
            self.change_commandline("c",tmptext)
        
        elif ( self.sender() == self.group_advanced.findChild(QtWidgets.QLineEdit,"line_s")):
            tmptext = self.group_advanced.findChild(QtWidgets.QLineEdit,"line_s").text()
            self.change_commandline("s",tmptext)
                
        elif ( self.sender() == self.group_advanced.findChild(QtWidgets.QLineEdit,"line_w")):
            tmptext = self.group_advanced.findChild(QtWidgets.QLineEdit,"line_w").text()
            self.change_commandline("w",tmptext)
        
        elif ( self.sender() == self.group_coupled.findChild(QtWidgets.QLineEdit,"line_PE_atmos") or
               self.sender() == self.group_coupled.findChild(QtWidgets.QLineEdit,"line_PE_ocean")):
            t1 = self.group_coupled.findChild(QtWidgets.QLineEdit,"line_PE_atmos").text()
            t2 = self.group_coupled.findChild(QtWidgets.QLineEdit,"line_PE_ocean").text()
            tmptext = "%s %s" % (t1,t2)
            
            self.change_commandline("Y",tmptext)

            # try catch sowas in der art
            pe1 = 0
            pe2 = 0
            
            try:
                pe1 = int(t1)
            except ValueError:            
                pass

            try:
                pe2 = int(t2)
            except ValueError:
                pass

            PE_total = pe1+pe2   
            self.group_execution.findChild(QtWidgets.QLineEdit,"line_pe").setText(str(PE_total))
            self.change_commandline("X",str(PE_total))

    # deletes parameter from commandline
    ##################################################################################### 
    def delete_commandline(self, id_str):   
        # Read palmrunline
        commandline = self.groupBox.findChild(QtWidgets.QLineEdit,"commandline")
        palmrunline = str(commandline.text())
        s = " -%s" % (id_str)
        splitline = list(filter(None,palmrunline.split(s)))
        
        if ( len(splitline) == 1):
            return 1
        else:
            param = splitline[1].split("-")
            param[0] = ""
            splitline[1]= " -".join(param)
            newpalmrunline = "".join(splitline)
            newpalmrunline = newpalmrunline.replace("  "," ")
    
            # Print new palmrunline to screen
            commandline.setText(newpalmrunline)
        
    # controls flags
    ###################################################################################
    def check_flags(self):
        status = self.group_execution.findChild(QtWidgets.QCheckBox,"check_delete_tmp_files" ).checkState()
        if (status == 2):      
            self.activate_flag("B")
        else:
            self.deactivate_flag("B")
    
        status = self.groupBox.findChild(QtWidgets.QCheckBox,"check_verbose" ).checkState()
        if (status == 2):
            self.activate_flag("v") 
        else:
            self.deactivate_flag("v")
        
        status = self.group_advanced.findChild(QtWidgets.QCheckBox,"check_b" ).checkState()
        if (status == 2):
            self.activate_flag("b")
        else:
            self.deactivate_flag("b")
    
        status = self.group_advanced.findChild(QtWidgets.QCheckBox,"check_F" ).checkState()
        if (status == 2):
            self.activate_flag("F")
        else:
            self.deactivate_flag("F")
            
        status = self.group_advanced.findChild(QtWidgets.QCheckBox,"check_I" ).checkState()
        if (status == 2):
            self.activate_flag("I")
        else:
            self.deactivate_flag("I")

        status = self.group_advanced.findChild(QtWidgets.QCheckBox,"check_k" ).checkState()
        if (status == 2):
            self.activate_flag("k")
        else:
            self.deactivate_flag("k")

        status = self.group_advanced.findChild(QtWidgets.QCheckBox,"check_O" ).checkState()
        if (status == 2):
            self.activate_flag("O")
        else:
            self.deactivate_flag("O")

        status = self.group_execution.findChild(QtWidgets.QCheckBox,"check_prebuilt" ).checkState()
        if (status == 2):
            self.activate_flag("V")
        else:
            self.deactivate_flag("V")

        status = self.group_advanced.findChild(QtWidgets.QCheckBox,"check_x" ).checkState()
        if (status == 2):
            self.activate_flag("x")
        else:
            self.deactivate_flag("x")
        
        status = self.group_advanced.findChild(QtWidgets.QCheckBox,"check_Z" ).checkState()
        if (status == 2):
            self.activate_flag("Z")
        else:
            self.deactivate_flag("Z")
        
    # changes flag to parameter 
    ##################################################################################    
    def activate_flag(self, id_str):
        commandline = self.groupBox.findChild(QtWidgets.QLineEdit,"commandline")
        palmrunline = str(commandline.text())
        s = " -%s" % (id_str)
        splitline = list(filter(None,palmrunline.split(s)))
        
        if ( len(splitline) == 1):
            splitline.append("")
            s = " -%s" % (id_str)
            newpalmrunline = s.join(splitline)
            newpalmrunline = newpalmrunline.replace("  "," ")
            commandline.setText(newpalmrunline)                
        
    # deletes flag in commandline
    ####################################################################################    
    def deactivate_flag(self, id_str):
        commandline = self.groupBox.findChild(QtWidgets.QLineEdit,"commandline")
        palmrunline = str(commandline.text())
        s = " -%s" % (id_str)
        splitline = list(filter(None,palmrunline.split(s)))
 
        newpalmrunline = "".join(splitline)
        newpalmrunline = newpalmrunline.replace("  "," ")
        commandline.setText(newpalmrunline)     

    # Clears window 
    #################################################################################
    def reset_window(self):
        self.setupUi(self)
        self.tabWidget.setCurrentIndex(0)
        self.palm_logo.setPixmap(QtGui.QPixmap(palm_dir + "/trunk/SCRIPTS/palmrungui_files/logo.png"))
        self.recent_jobs(50)
        self.load_jobs()
        

    # Safes current commandline and user parameters to default file
    ################################################################################
    def save_default(self):
        string_to_file = self.groupBox.findChild(QtWidgets.QLineEdit,"commandline").text()

        filename ="%s/.palmrungui.default" % (palm_dir)
        tmstamp = strftime("%Y/%m/%d %H:%M")

        file = open(filename,"w")
        s = "%s %s" % (tmstamp, string_to_file)
        file.write(s)
        file.close()        

        
    # Executes command which starts watchdog (start palm_wd)
    ########################################################
    def start_watchdog(self):
        subprocess.Popen(["nohup","palm_wd",">>","/dev/null", "2>&1","&"])

    # Opens "help" dialog
    #################################
    def help(self):        
        dialog = HelpDialog()
        dialog.exec_()

    # Opens "about" dialog
    ##################################
    def about_gui(self):
        dialog = AboutDialog()
        dialog.exec_()

    # commandline to buttons etc
    ##############################
    def setup_gui(self, palmrun_str):  
      
        #self.palm_logo.setPixmap(QtWidgets.QPixmap(palm_dir + "/trunk/SCRIPTS/palmrungui_files/logo.png"))
      
        #  Some initial settings
        coupled_run = False 
        ocean_run   = False 
        nojob       = False 

        #Split parameters in palmrunline
        splitline = palmrun_str.split(" -")
        
        if ( splitline[0] != "palmrun"):
            return 1

        else:
            self.group_execution.setEnabled(True) 

        # Loop for all parameters in palmrunline
        i = len(splitline)-1
        while i >= 1:

            # Determine parameter
            splitparameter = splitline[i].split(" ")

            parameter = splitparameter[0] 
            splitparameter.pop(0)
            options = " ".join(splitparameter) 
            options = options.replace("\"","") 

            # Check for suitable switch
            if ( parameter == "r"):
              
                if ( options != ""):
                    self.group_execution.findChild(QtWidgets.QLineEdit,"line_jobname").setText(options)
                    nojob = False 
                
                else:                 
                    nojob = True 
                
            elif ( parameter == "c"):
                tmpindex = self.group_execution.findChild(QtWidgets.QComboBox, "combo_configuration").findText(options.strip(),QtCore.Qt.MatchExactly)
                if tmpindex != -1: 
                   self.group_execution.findChild(QtWidgets.QComboBox, "combo_configuration").setCurrentIndex(tmpindex)
                else:
                   self.group_execution.findChild(QtWidgets.QComboBox, "combo_configuration").setCurrentIndex(0)
            elif ( parameter == "q"):
                self.group_execution.findChild(QtWidgets.QLineEdit,"line_q").setText(options)
             
            elif ( parameter == "A"):
                self.group_execution.findChild(QtWidgets.QLineEdit,"line_account").setText(options)

            elif ( parameter == "W"):
                self.group_execution.findChild(QtWidgets.QLineEdit,"line_wait").setText(options)

            elif ( parameter == "X"):
                self.group_execution.findChild(QtWidgets.QLineEdit,"line_pe").setText(options)
             
            elif ( parameter == "T"):
                self.group_execution.findChild(QtWidgets.QLineEdit,"line_tpn").setText(options)
             
            elif ( parameter == "t"):
                self.group_execution.findChild(QtWidgets.QLineEdit,"line_time").setText(options)
             
            elif ( parameter == "B"):
                self.group_execution.findChild(QtWidgets.QCheckBox,"check_delete_tmp_files").setChecked(True)
             
            elif ( parameter == "v"):
                self.groupBox.findChild(QtWidgets.QCheckBox,"check_verbose").setChecked(True)
                          
            elif ( parameter == "b"): 
                self.group_advanced.findChild(QtWidgets.QCheckBox,"check_b").setChecked(True)
             
            elif ( parameter == "F"):
                self.group_advanced.findChild(QtWidgets.QCheckBox,"check_F").setChecked(True)
             
            elif ( parameter == "I"):
                self.group_advanced.findChild(QtWidgets.QCheckBox,"check_I").setChecked(True)
             
            elif ( parameter == "k"):
                self.group_advanced.findChild(QtWidgets.QCheckBox,"check_k").setChecked(True)
             
            elif ( parameter == "O"): 
                self.group_advanced.findChild(QtWidgets.QCheckBox,"check_O").setChecked(True)

            elif ( parameter == "V"):
                self.group_execution.findChild(QtWidgets.QCheckBox,"check_prebuilt").setChecked(True)

            elif ( parameter == "x"):              
                self.group_advanced.findChild(QtWidgets.QCheckBox,"check_x").setChecked(True)

            elif ( parameter == "Z"):
                self.group_advanced.findChild(QtWidgets.QCheckBox,"check_Z").setChecked(True)
               
            elif ( parameter == "m"):
                self.group_advanced.findChild(QtWidgets.QLineEdit,"line_m").setText(options)
             
            elif ( parameter == "M"):
                self.group_advanced.findChild(QtWidgets.QLineEdit,"line_M").setText(options)
                          
            elif ( parameter == "D"):
                self.group_advanced.findChild(QtWidgets.QLineEdit,"line_D").setText(options)
             
            elif ( parameter == "c"):
                self.group_advanced.findChild(QtWidgets.QLineEdit,"line_c").setText(options)
                          
            elif ( parameter == "s"):
                self.group_advanced.findChild(QtWidgets.QLineEdit,"line_s").setText(options)
             
            elif ( parameter == "w"):
                self.group_advanced.findChild(QtWidgets.QLineEdit,"line_w").setText(options)
             

            # Determine settings for coupled restart runs
            elif ( parameter == "Y"):
                optionssplit = options.split(" ") 
                
                group_coupled = self.group_execution.findChild(QtWidgets.QGroupBox,"group_coupled")
                group_coupled.findChild(QtWidgets.QLineEdit,"line_PE_atmos").setEnabled(True)
                group_coupled.findChild(QtWidgets.QLineEdit,"line_PE_ocean").setEnabled(True)
                group_coupled.findChild(QtWidgets.QLabel,"label_coupled1").setEnabled(True)
                group_coupled.findChild(QtWidgets.QLabel,"label_coupled2").setEnabled(True)
                group_coupled.findChild(QtWidgets.QLabel,"label_coupled3").setEnabled(True)
                group_coupled.findChild(QtWidgets.QLabel,"label_coupling").setEnabled(True)

                if (optionssplit.count() == 2):
                    group_coupled.findChild(QtWidgets.QLineEdit,"line_PE_atmos").setEnabled(optionssplit[0])
                    group_coupled.findChild(QtWidgets.QLineEdit,"line_PE_ocean").setEnabled(optionssplit[1])
                 
                else:                  
                    group_coupled.findChild(QtWidgets.QLineEdit,"line_PE_atmos").setText("")
                    group_coupled.findChild(QtWidgets.QLineEdit,"line_PE_ocean").setText("")
 
                 
                coupled_run = True 
            
            elif ( parameter == "y"):              
                self.group_execution.findChild(QtWidgets.QComboBox, "drop_job").setCurrentIndex(3)
             

            # Determine settings for the run control list
            elif ( parameter == "a"):
              
                optionssplit = options.split(" ") 

                options_2 = None
                options_all = None
                
                j = 0
                while j < len(optionssplit):
                 
                    options_all = optionssplit[j] 
                    options_2 = optionssplit[j][:2] 
                                              
                    if (options_2 == "d3"):      
                        if (options_all[:3][-1] == "#"):
                            self.group_execution.findChild(QtWidgets.QComboBox, "drop_job").setCurrentIndex(0)
                        elif (options_all[:3][-1] == "r"):
                            self.group_execution.findChild(QtWidgets.QComboBox, "drop_job").setCurrentIndex(1)
                        
                        elif (options_all[:3][-1] == "o"):
                            ocean_run = True 
                        
                    if (options_all == "restart"):
                        self.group_execution.findChild(QtWidgets.QCheckBox,"check_restarts").setChecked(True)

                    if (options_all == "rec"):
                        self.group_execution.findChild(QtWidgets.QCheckBox,"check_cycfill").setChecked(True)

                    if (options_all == "svfin"):
                        self.group_execution.findChild(QtWidgets.QCheckBox,"check_svf").setChecked(True)

                    if (options_all == "svfout"):
                        self.group_execution.findChild(QtWidgets.QCheckBox,"check_svf_2").setChecked(True)

                    j = j+1
            i = i-1
        # Change drop box state in case of ocean precursor or coupled restart runs
        if ( ocean_run == True ):
            if ( coupled_run == True ):
                self.group_execution.findChild(QtWidgets.QComboBox, "drop_job").setCurrentIndex(4)
            
            else:
                self.group_execution.findChild(QtWidgets.QComboBox, "drop_job").setCurrentIndex(3)

        # Join palmrunline and post it to mainwindow
        palmrunline = " -".join(splitline) 
        self.groupBox.findChild(QtWidgets.QLineEdit,"commandline").setText(palmrunline)

        # Disable mainwindow if no job was found, otherwise enable
        if ( nojob == True ):
            self.group_execution.setEnabled(False) 
            self.groupBox.findChild(QtWidgets.QPushButton,"button_start").setEnabled(False)
            self.menuBar.findChild(QtWidgets.QMenu,"menuStart").actions()[3].setEnabled(False)
            self.group_execution.findChild(QtWidgets.QComboBox, "drop_job").setEnabled(False)
            self.group_advanced.setEnabled(False) 
            self.check_advanced.setEnabled(False) 
            self.groupBox.findChild(QtWidgets.QLineEdit,"line_tag").setEnabled(False)
         
        else:
            self.group_execution.setEnabled(True) 
            self.groupBox.findChild(QtWidgets.QPushButton,"button_start").setEnabled(True)
            self.menuBar.findChild(QtWidgets.QMenu,"menuStart").actions()[3].setEnabled(True)
            self.group_execution.findChild(QtWidgets.QComboBox, "drop_job").setEnabled(True)
            self.group_advanced.setEnabled(True) 
            self.groupBox.findChild(QtWidgets.QLineEdit,"line_tag").setEnabled(True)


        self.tabWidget.setCurrentIndex(0)


    ## Open from history
    ###################################
    def open_from_file(self):

#      Show History
       opt = OpenHistoryBox()
       opt.exec_()

       if ( history_entry != "" ):
          palmrunline = str(history_entry)
          palmrunline = palmrunline[17:]
          palmrunline = palmrunline[:len(palmrunline)-1]
          palmrunline = list(filter(None,palmrunline.split("(")))[0]
          self.groupBox.findChild(QtWidgets.QLineEdit,"commandline").setText(palmrunline)

#         Set selected item to jobname in "available jobs" list
          jobname = str(history_entry[17:])
          jobname = list(filter(None,jobname.split(" -r")))[1]
          jobname = jobname.strip()
          jobname = list(filter(None,jobname.split(" ")))[0]
          jobname = jobname.replace("\"","")
          item2int = self.list_jobs.findItems(jobname,QtCore.Qt.MatchCaseSensitive)

          if ( item2int != [] ):
             self.list_jobs.setCurrentItem(item2int[0])
             self.update_input()

#         Add tooltip tag
          tag = str(history_entry).split('\n')[0]
          tag = list(filter(None,tag.split("(")))[1]
          tag = tag.replace(")","")
          self.groupBox.findChild(QtWidgets.QLineEdit,"line_tag").setText(tag)
 
#         Process palmrungui to set up gui controls
          self.setup_gui(palmrunline)


    def update_all(self):
 
       self.setEnabled(False)
       self.list_input.clear() 
       self.list_user.clear()
       self.list_monitoring.clear()
       self.list_output.clear()
       self.load_jobs()
       self.setEnabled(True)
       self.update_input()


    # Load jobs into list
    def load_jobs(self):
       
       selected_job = self.list_jobs.currentItem()
       
       if ( selected_job is not None ):
          jobname = selected_job.text()
       else:
          jobname = ""

       self.list_jobs.clear() 
       self.line_path.setText(job_dir + "/")
       
       list_of_files = os.listdir(job_dir)
      

       for i in range(0,len(list_of_files)):
          tmp_file = job_dir + "/" + list_of_files[i]

          if ( os.path.isdir(tmp_file) ):
             self.list_jobs.addItem(str(list_of_files[i]))
       


       item2int = self.list_jobs.findItems(jobname,QtCore.Qt.MatchCaseSensitive)
        
       if ( item2int != [] ):
          self.list_jobs.setCurrentItem(item2int[0])


    # Update input and user code lists
    def update_input(self):

       self.labeltimer.stop() 

      
       self.list_input.clear() 
       self.list_user.clear()
       self.list_monitoring.clear()
       self.list_output.clear()
       
       jobitem = self.list_jobs.currentItem()
       
       if ( jobitem != None ):
       
          job_to_show = job_dir + "/" + jobitem.text() + "/INPUT"

          if ( os.path.isdir(job_to_show) ):

             list_of_files = os.listdir(job_to_show)
         
             for i in range(0,len(list_of_files)):
                tmp_file = job_to_show + "/" + list_of_files[i]
             
                if ( os.path.isfile(tmp_file) ):
                   self.list_input.addItem(str(list_of_files[i]) + " (" + file_size(tmp_file) + ")")
  
          job_to_show = job_dir + "/" + jobitem.text() + "/USER_CODE"
       
          if ( os.path.isdir(job_to_show) ):
          
             list_of_files = os.listdir(job_to_show)
          
             for i in range(0,len(list_of_files)):
                tmp_file = job_to_show + "/" + list_of_files[i]

                if ( os.path.isfile(tmp_file) ):
                   self.list_user.addItem(str(list_of_files[i]))

          job_to_show = job_dir + "/" + jobitem.text() + "/MONITORING"
       
          if ( os.path.isdir(job_to_show) ):
          
             list_of_files = os.listdir(job_to_show)
          
             for i in range(0,len(list_of_files)):
                tmp_file = job_to_show + "/" + list_of_files[i]

                if ( os.path.isfile(tmp_file) ):
                   self.list_monitoring.addItem(str(list_of_files[i]) + " (" + file_size(tmp_file) + ")")

          job_to_show = job_dir + "/" + jobitem.text() + "/OUTPUT"

          if ( os.path.isdir(job_to_show) ):

             list_of_files = os.listdir(job_to_show)

             for i in range(0,len(list_of_files)):
                tmp_file = job_to_show + "/" + list_of_files[i]

                if ( os.path.isfile(tmp_file) ):
                   self.list_output.addItem(str(list_of_files[i]) + " (" + file_size(tmp_file) + ")")
  
          self.group_execution.findChild(QtWidgets.QLineEdit,"line_jobname").setText(jobitem.text())
       self.group_history.findChild(QtWidgets.QListWidget,"list_jobname").clearSelection()

       self.timetimer.start()
       self.timer.start(update_frequency)
       self.labeltimer.start(int(update_frequency/10))
       self.push_update.setText("Update (" + str(int(update_frequency/1000/60)) + " min)")
       QtWidgets.QApplication.processEvents()

       # Change palmrunline accordingly
       self.change_commandline("r","")
       self.change_commandline("a","")  

    # Create a new (empty) job
    def create_new_job(self):

       self.setEnabled(False)

       text, ret = QtWidgets.QInputDialog.getText(self, "Create new job", "Enter new job name:", QtWidgets.QLineEdit.Normal, text = '')

       if ( ret ):
          new_job_name = str(text)
       else:
          self.setEnabled(True)
          return

       new_input_dir  = job_dir + "/" + new_job_name + "/INPUT"

#      check if a job exists with the new job name. If not, the new job is created and an empty _p3d file is created.
       if ( os.path.isdir(new_input_dir) ):
          notify = QtWidgets.QMessageBox.warning(self,'Create new job directory',"Error. Could not create job directory. A job with the new name already exists.")
          self.setEnabled(True)
          return
       else:
          os.makedirs(new_input_dir)
          os.mknod(new_input_dir + '/' + new_job_name + '_p3d')

#      write default namelists into _p3d file
       filename = open(new_input_dir + '/' + new_job_name + '_p3d', "w")
       filename.write('&initialization_parameters\n    nx = 39,\n    ny = 39,\n    nz = 40,\n    dx = 50.0,\n    dy = 50.0,\n    dz = 50.0,\n/\n\n')
       filename.write('&runtime_parameters\n    end_time = 3600.0,\n/\n')
       filename.close()


#      create empty USER_CODE folder for convenience
       new_user_dir  = job_dir + "/" + new_job_name + "/USER_CODE"
       os.makedirs(new_user_dir)


       self.load_jobs()
       self.list_input.clear()
       self.list_user.clear()
       self.list_monitoring.clear()
       self.list_output.clear()
       self.setEnabled(True)


    # Make a copy of a job
    def copy_job(self):

       self.setEnabled(False)
       old_job_name = self.list_jobs.currentItem().text()
       
       text, ret = QtWidgets.QInputDialog.getText(self, "Copy job", "Enter new job name:", QtWidgets.QLineEdit.Normal, text = old_job_name)
       if ( ret ):
          new_job_name = str(text)
       else:
          self.setEnabled(True)
          return

       new_input_dir  = job_dir + "/" + new_job_name + "/INPUT"

#      check if a job exists with the new job name  
       if ( os.path.isdir(new_input_dir) ):
          notify = QtWidgets.QMessageBox.warning(self,'Create new job directory',"Error. Could not create job directory. A job with the new name already exists.")
          self.setEnabled(True)
          return
       else:
          os.makedirs(new_input_dir)

#      copy and rename input files (if present)       
       job_to_copy = job_dir + "/" + old_job_name + "/INPUT"

       if ( os.path.isdir(job_to_copy) ):

          list_of_files = os.listdir(job_to_copy)

          for i in range(0,len(list_of_files)):

             tmp_file = job_to_copy + "/" + list_of_files[i]          
             new_file = new_input_dir + "/" + list_of_files[i].replace(old_job_name, new_job_name)
             shutil.copy(tmp_file, new_file)




       new_user_dir  = job_dir + "/" + new_job_name + "/USER_CODE"

#      check if user code exists in the new job directory 
       if ( os.path.isdir(new_user_dir) ):
          notify = QtWidgets.QMessageBox.warning(self,'Create new user code directory',"Error. Could not create user code directory. A user code directiory with the new name already exists.")
          self.setEnabled(True)
          return
       else:
          os.makedirs(new_user_dir)


#      copy user code files (if present)       
       user_to_copy = job_dir + "/" + old_job_name + "/USER_CODE"

       if ( os.path.isdir(user_to_copy) ):

          list_of_files = os.listdir(user_to_copy)
         
          for i in range(0,len(list_of_files)):

             tmp_file = user_to_copy + "/" + list_of_files[i]          
             new_file = new_user_dir + "/" + list_of_files[i]
             shutil.copy(tmp_file, new_file)

       self.load_jobs()
       self.list_input.clear() 
       self.list_user.clear()
       self.list_monitoring.clear()
       self.list_output.clear()
       self.setEnabled(True)


    # Create a whole set of jobs
    def create_set(self):
  
       global set_list
#      disable mainwindow  
       self.setEnabled(False)
      
#      show Options Dialog      
       opt = CreateSetBox()
       opt.exec_()

       old_job_name = self.list_jobs.currentItem().text()

       for j in range(0,len(set_list)):

          if ( set_list[j] != "" ):
             new_job_name   = str(set_list[j])
             new_input_dir  = job_dir + "/" + str(set_list[j]) + "/INPUT"
          else:
             continue

#         check if a job exists with the new job name  
          if ( os.path.isdir(new_input_dir) ):
             notify = QtWidgets.QMessageBox.warning(self,'Create new job directory',"Error. Could not create job directory. A job with the new name already exists.")
             self.setEnabled(True)
             return
          else:
             os.makedirs(new_input_dir)

#         copy and rename input files (if present)       
          job_to_copy = job_dir + "/" + old_job_name + "/INPUT"

          if ( os.path.isdir(job_to_copy) ):

             list_of_files = os.listdir(job_to_copy)
         
             for i in range(0,len(list_of_files)):

                tmp_file = job_to_copy + "/" + list_of_files[i]          
                new_file = new_input_dir + "/" + list_of_files[i].replace(old_job_name, new_job_name)
                shutil.copy(tmp_file, new_file)


          new_user_dir  = job_dir + "/" + new_job_name + "/USER_CODE"

#         check if user code exists in the new job directory 
          if ( os.path.isdir(new_user_dir) ):
             notify = QtWidgets.QMessageBox.warning(self,'Create new user code directory',"Error. Could not create user code directory. A user code directiory with the new name already exists.")
             self.setEnabled(True)
             return
          else:
             os.makedirs(new_user_dir)


#         copy user code files (if present)       
          user_to_copy = job_dir + "/" + old_job_name + "/USER_CODE"

          if ( os.path.isdir(user_to_copy) ):

             list_of_files = os.listdir(user_to_copy)
         
             for i in range(0,len(list_of_files)):

                tmp_file = user_to_copy + "/" + list_of_files[i]          
                new_file = new_user_dir + "/" + list_of_files[i]
                shutil.copy(tmp_file, new_file)

          self.load_jobs()
          self.list_input.clear() 
          self.list_user.clear()
          self.list_monitoring.clear()
          self.list_output.clear()

       self.setEnabled(True) 
       set_list = []


 # Add a custom context menu for the job selection list
    def openmenujob(self, position):

        menu = QtWidgets.QMenu()


        newAction = QtWidgets.QAction('Create new job', self)
        newAction.triggered.connect(self.create_new_job)
        menu.addAction(newAction)
        selection = self.list_jobs.selectedItems()

        if ( len(selection) != 0 ):


           copyAction = QtWidgets.QAction('Copy job', self)
           copyAction.triggered.connect(self.copy_job)
           createAction = QtWidgets.QAction('Create set from job', self)
           createAction.triggered.connect(self.create_set)
           delAction = QtWidgets.QAction('Delete job', self)
           delAction.triggered.connect(self.DeleteJob)

           menu.addAction(copyAction)
           menu.addAction(createAction)
           menu.addAction(delAction)

        action = menu.exec_(self.list_jobs.mapToGlobal(position))

 # Add a custom context menu
    def openmenuinput(self, position):

        menu = QtWidgets.QMenu()

        selection = self.list_input.selectedItems()
        
        if ( len(selection) != 0 ):


           openAction = QtWidgets.QAction('Open file(s)', self)
           openAction.setStatusTip('Open file(s) in your favorite editor')
           openAction.triggered.connect(self.OpenFilesInput)


           delAction = QtWidgets.QAction('Delete selected file(s)', self)
           delAction.triggered.connect(self.DeleteFilesInput)

           menu.addAction(openAction)
           menu.addAction(delAction)

        delallAction = QtWidgets.QAction('Delete all files', self)
        delallAction.triggered.connect(self.DeleteAllFilesInput)
        menu.addAction(delallAction)

        action = menu.exec_(self.list_input.mapToGlobal(position))

 # Add a custom context menu
    def openmenuuser(self, position):

        menu = QtWidgets.QMenu()

        selection = self.list_user.selectedItems()
        
        if ( len(selection) != 0 ):

           openAction = QtWidgets.QAction('Open file(s)', self)
           openAction.setStatusTip('Open file(s) in your favorite editor')
           openAction.triggered.connect(self.OpenFilesUser)


           delAction = QtWidgets.QAction('Delete file(s)', self)
           delAction.triggered.connect(self.DeleteFilesUser)

           menu.addAction(openAction)
           menu.addAction(delAction)

        delallAction = QtWidgets.QAction('Delete all files', self)
        delallAction.triggered.connect(self.DeleteAllFilesUser)
        menu.addAction(delallAction)

        action = menu.exec_(self.list_user.mapToGlobal(position))

 # Add a custom context menu
    def openmenuoutput(self, position):

        menu = QtWidgets.QMenu()

        selection = self.list_output.selectedItems()

        if ( len(selection) != 0 ):

           openAction = QtWidgets.QAction('Open file(s)', self)
           openAction.setStatusTip('Open file(s) in your favorite editor')
           openAction.triggered.connect(self.OpenFilesOutput)

           delAction = QtWidgets.QAction('Delete file(s)', self)
           delAction.triggered.connect(self.DeleteFilesOutput)



           menu.addAction(openAction)
           menu.addAction(delAction)

        delallAction = QtWidgets.QAction('Delete all files', self)
        delallAction.triggered.connect(self.DeleteAllFilesOutput)
        menu.addAction(delallAction)

        action = menu.exec_(self.list_output.mapToGlobal(position))

 # Add a custom context menu
    def openmenumonitoring(self, position):

        menu = QtWidgets.QMenu()

        selection = self.list_monitoring.selectedItems()
        
        if ( len(selection) != 0 ):

           openAction = QtWidgets.QAction('Open file(s)', self)
           openAction.setStatusTip('Open file(s) in your favorite editor')
           openAction.triggered.connect(self.OpenFilesMonitoring)

           delAction = QtWidgets.QAction('Delete file(s)', self)
           delAction.triggered.connect(self.DeleteFilesMonitoring)



           menu.addAction(openAction)
           menu.addAction(delAction)

        delallAction = QtWidgets.QAction('Delete all files', self)
        delallAction.triggered.connect(self.DeleteAllFilesMonitoring)
        menu.addAction(delallAction)

        action = menu.exec_(self.list_monitoring.mapToGlobal(position))

    def OpenConfig(self):

       config = str(self.group_execution.findChild(QtWidgets.QComboBox, "combo_configuration").currentText())
       if ( config != "" ):
          filename = ".palm.config." + config
          open_file = "xdg-open " + filename
          os.system(str(open_file))

    def OpenFilesInput(self):
   
       sel_job = self.list_jobs.currentItem().text()
       sel_files = self.list_input.selectedItems()
       
       input_dir = job_dir + "/" + sel_job + "/INPUT/"
       
       open_files = ""
       for i in range(0,len(sel_files)):
          open_files = open_files + "xdg-open " + input_dir +  sel_files[i].text().split("(")[0] + "; "

       os.system(str(open_files))

    def OpenFilesUser(self):
   
       sel_job = self.list_jobs.currentItem().text()
       sel_files = self.list_user.selectedItems()
       
       input_dir = job_dir + "/" + sel_job + "/USER_CODE/"
       
       open_files = ""
       for i in range(0,len(sel_files)):
          open_files = open_files + "xdg-open " + input_dir +  sel_files[i].text() + "; "

       os.system(str(open_files)) 
       selection = self.list_jobs.selectedItems()
  
    def OpenFilesMonitoring(self):
   
       sel_job = self.list_jobs.currentItem().text()
       sel_files = self.list_monitoring.selectedItems()
       
       input_dir = job_dir + "/" + sel_job + "/MONITORING/"
       
       open_files = ""
       for i in range(0,len(sel_files)):
          open_files = open_files + "xdg-open " + input_dir +  sel_files[i].text().split("(")[0]  + "; "

       os.system(str(open_files)) 
       selection = self.list_jobs.selectedItems()
  
    def OpenFilesOutput(self):
   
       sel_job = self.list_jobs.currentItem().text()
       sel_files = self.list_output.selectedItems()
       
       input_dir = job_dir + "/" + sel_job + "/OUTPUT/"
       
       open_files = ""
       for i in range(0,len(sel_files)):
          open_files = open_files + "xdg-open " + input_dir +  sel_files[i].text().split("(")[0] + "; "

       os.system(str(open_files)) 
       selection = self.list_jobs.selectedItems()

    def DeleteFilesInput(self):

       status = QtWidgets.QMessageBox.question(self,'Delete Confirmation', "Are you sure you want to delete these files?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

       if status == QtWidgets.QMessageBox.Yes:

          sel_job = self.list_jobs.currentItem().text()
          sel_files = self.list_input.selectedItems()

          input_dir = job_dir + "/" + sel_job + "/INPUT/"

          for i in range(0,len(sel_files)):
             filename = input_dir +  sel_files[i].text().split("(")[0].strip()
             os.remove(filename)

          self.update_all()


    def DeleteFilesUser(self):

       status = QtWidgets.QMessageBox.question(self,'Delete Confirmation', "Are you sure you want to delete these files?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

       if status == QtWidgets.QMessageBox.Yes:

          sel_job = self.list_jobs.currentItem().text()
          sel_files = self.list_user.selectedItems()

          input_dir = job_dir + "/" + sel_job + "/USER_CODE/"

          for i in range(0,len(sel_files)):
             filename = input_dir +  sel_files[i].text().split("(")[0].strip()
             os.remove(filename)

          self.update_all()


    def DeleteFilesMonitoring(self):

       status = QtWidgets.QMessageBox.question(self,'Delete Confirmation', "Are you sure you want to delete these files?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

       if status == QtWidgets.QMessageBox.Yes:

          sel_job = self.list_jobs.currentItem().text()
          sel_files = self.list_monitoring.selectedItems()

          input_dir = job_dir + "/" + sel_job + "/MONITORING/"

          for i in range(0,len(sel_files)):
             filename = input_dir +  sel_files[i].text().split("(")[0].strip()
             os.remove(filename)

          self.update_all()


    def DeleteFilesOutput(self):

       status = QtWidgets.QMessageBox.question(self,'Delete Confirmation', "Are you sure you want to delete these files?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

       if status == QtWidgets.QMessageBox.Yes:

          sel_job = self.list_jobs.currentItem().text()
          sel_files = self.list_output.selectedItems()

          input_dir = job_dir + "/" + sel_job + "/OUTPUT/"

          for i in range(0,len(sel_files)):
             filename = input_dir +  sel_files[i].text().split("(")[0].strip()
             os.remove(filename)

          self.update_all()


    def DeleteAllFilesInput(self):

       status = QtWidgets.QMessageBox.question(self,'Delete Confirmation', "Are you sure you want to delete all files?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

       if status == QtWidgets.QMessageBox.Yes:

          check = self.list_jobs.currentItem()
          if ( check ):
             sel_job = self.list_jobs.currentItem().text()
             sel_files = self.list_input.selectedItems()

             input_dir = job_dir + "/" + sel_job + "/INPUT/"

             if os.path.exists(str(input_dir)):
                shutil.rmtree(str(input_dir))
                os.makedirs(str(input_dir))

                self.update_all()


    def DeleteAllFilesUser(self):

       status = QtWidgets.QMessageBox.question(self,'Delete Confirmation', "Are you sure you want to delete all files?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

       if status == QtWidgets.QMessageBox.Yes:

          check = self.list_jobs.currentItem()
          if ( check ):
             sel_job = self.list_jobs.currentItem().text()
             sel_files = self.list_input.selectedItems()

             input_dir = job_dir + "/" + sel_job + "/USER_CODE/"

             if os.path.exists(str(input_dir)):
                shutil.rmtree(str(input_dir))
                os.makedirs(str(input_dir))

                self.update_all()



    def DeleteAllFilesOutput(self):

       status = QtWidgets.QMessageBox.question(self,'Delete Confirmation', "Are you sure you want to delete all files?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

       if status == QtWidgets.QMessageBox.Yes:

          check = self.list_jobs.currentItem()
          if ( check ):
             sel_job = self.list_jobs.currentItem().text()
             sel_files = self.list_input.selectedItems()

             input_dir = job_dir + "/" + sel_job + "/OUTPUT/"

             if os.path.exists(str(input_dir)):
                shutil.rmtree(str(input_dir))
                os.makedirs(str(input_dir))

                self.update_all()



    def DeleteAllFilesMonitoring(self):

       status = QtWidgets.QMessageBox.question(self,'Delete Confirmation', "Are you sure you want to delete all files?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

       if status == QtWidgets.QMessageBox.Yes:

          check = self.list_jobs.currentItem()
          if ( check ):
             sel_job = self.list_jobs.currentItem().text()
             sel_files = self.list_input.selectedItems()

             input_dir = job_dir + "/" + sel_job + "/MONITORING/"

             if os.path.exists(str(input_dir)):
                shutil.rmtree(str(input_dir))
                os.makedirs(str(input_dir))

                self.update_all()



    def DeleteJob(self):

       status = QtWidgets.QMessageBox.question(self,'Delete Confirmation', "Are you sure you want to delete this job?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

       if status == QtWidgets.QMessageBox.Yes:

          sel_job = self.list_jobs.currentItem().text()

          input_dir = job_dir + "/" + sel_job
          shutil.rmtree(str(input_dir))

          self.update_all()

# Message box for showing RUN_CONTROL output
class CreateSetBox(QtWidgets.QDialog):
    def __init__(self):

        super(CreateSetBox, self).__init__()

        uic.loadUi(palmrungui_sdir + '/palmrungui_files/create_set.ui', self)

        self.show()

        return

#   Cancel button
    def rejected(self):

       self.close()

       return

#   OK button
    def accept(self):

       global set_list

       text = self.list.toPlainText() 
       set_list = text.split('\n')
       self.close()

       return

# Message box for showing RUN_CONTROL output
class OpenHistoryBox(QtWidgets.QDialog):
    def __init__(self):

        super(OpenHistoryBox, self).__init__()

        uic.loadUi(palmrungui_sdir + '/palmrungui_files/history.ui', self)

        if os.path.exists(palm_dir + "/.palmrungui.history"):
            pass
        else:
            return

        filename = open(palm_dir + "/.palmrungui.history","r")
        history = filename.readlines()
        filename.close()


        list_jobname = self.findChild(QtWidgets.QListWidget,"list_history")

        # Read history entries and append to recent job list
        len_history=len(history)-1
        i = 0
        while i<=len_history:
            list_jobname.addItem(history[i])
            i = i + 1

        self.show()

        return

#   Select item / activate Load button
    def ItemSelected(self):

       self.push_open.setEnabled(True)

       return


#   Load job from history
    def OpenFromHistory(self):

       global history_entry
       history_entry = self.list_history.selectedItems()[0].text()

       self.close()

       return

#   Close history window
    def DiscardHistory(self):

       global history_entry
       history_entry = ""
       self.close()

       return

#   Clear history
    def ClearHistory(self):

       status = QtWidgets.QMessageBox.question(self,'Delete Confirmation', "Are you sure you want to delete your history?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

       if status == QtWidgets.QMessageBox.Yes:

          if os.path.exists(palm_dir + "/.palmrungui.history"):
             os.remove(palm_dir + "/.palmrungui.history")


       return

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = Mainwindow()
    window.show()
    sys.exit(app.exec_())
