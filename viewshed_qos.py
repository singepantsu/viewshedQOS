# -*- coding: utf-8 -*-
"""
/***************************************************************************
 ViewshedQOS
                                 A QGIS plugin
 Run viewshed analysis based on our requirements.
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2023-06-26
        git sha              : $Format:%H$
        copyright            : (C) 2023 by Mekki Benhamidouche
        email                : mekki99@hotmail.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from qgis.PyQt.QtCore import QSettings, QTranslator, QCoreApplication
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction
from qgis.core import QgsProject
from .code.functionslist import *

# Initialize Qt resources from file resources.py
from .resources import *
# Import the code for the dialog
from .viewshed_qos_dialog import ViewshedQOSDialog
from .viewshed_qos_dialog_hardware import ViewshedQOSDialogHardware
import os.path
from osgeo import gdal


class ViewshedQOS:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'ViewshedQOS_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)
            QCoreApplication.installTranslator(self.translator)

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&Viewshed QOS')

        # Check if plugin was started the first time in current QGIS session
        # Must be set in initGui() to survive plugin reloads
        self.first_start = None

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('ViewshedQOS', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            # Adds plugin icon to Plugins toolbar
            self.iface.addToolBarIcon(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/viewshed_qos/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'Launch viewshed plugin'),
            callback=self.run,
            parent=self.iface.mainWindow())

        # will be set False in run()
        self.first_start = True


    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&Viewshed QOS'),
                action)
            self.iface.removeToolBarIcon(action)

    def open_hardware_windows(self):
        self.dlg.close()
        self.dlg = ViewshedQOSDialogHardware()
        
    def print_results(self):
        lumen = self.dlg.spinbox.value()
        candela = self.dlg.spinbox.value()
        
        # res = "Luminous flux:" + str(lumen) +"Lm" + '\n' + "Beam Power:" + str(candela) + 
        # self.dlg.textBrowser.setText(res)
    
    def addLumenValue(self):
        self.dlg.spinBox.setValue(calculateLumen(self.dlg.spinBox_3.value(),self.dlg.spinBox_2.value()))
    
    def addCandelaValue(self):
        self.dlg.spinBox_2.setValue(calculateCandela(self.dlg.spinBox.value(),self.dlg.spinBox_3.value()))
    
    def addBeamAngleValue(self):
        self.dlg.spinBox_3.setValue(calculateBeamAngle(self.dlg.spinBox.value(),self.dlg.spinBox_2.value()))
        
    def createRaster(self):
        self.dlg.runButton.setEnabled(False)
        selectedReliefLayer = (QgsProject.instance().mapLayersByName(self.dlg.comboBox.currentText())[0])
        px_size = selectedReliefLayer.rasterUnitsPerPixelX()
        selectedReliefLayer = selectedReliefLayer.dataProvider()
        vis_lay = (QgsProject.instance().mapLayersByName(self.dlg.comboBox_6.currentText())[0]).dataProvider()
        selectedSlopeLayer = QgsProject.instance().mapLayersByName(self.dlg.comboBox_5.currentText())[0]
        selectedBuildingLayer = (QgsProject.instance().mapLayersByName(self.dlg.comboBox_7.currentText())[0])
        selectedTarget = QgsProject.instance().mapLayersByName(self.dlg.comboBox_2.currentText())[0]
        selectedSource = QgsProject.instance().mapLayersByName(self.dlg.comboBox_3.currentText())[0]
        # sourceFeats = list(selectedSource.getFeatures())
        # target_coord = list(map(lambda x: (x.geometry().asPoint().x(),x.geometry().asPoint().y()),selectedTarget.getFeatures()))
        # source_coord = list(map(lambda x: (x.geometry().asPoint().x(),x.geometry().asPoint().y()),selectedSource.getFeatures()))
        target_coord = list(map(lambda x: x.geometry().asPoint(),selectedTarget.getFeatures()))
        source_coord = list(map(lambda x: x.geometry().asPoint(),selectedSource.getFeatures()))
        # source_coord = []
        # for src in sourceFeats:
        #     source_coord.append(src.geometry().asPoint())
        lumen = self.dlg.spinBox.value()
        beamAngle = self.dlg.spinBox_3.value()
        source_height = self.dlg.spinBox_4.value() + selectedReliefLayer.sample(QgsPointXY(source_coord[0][0],source_coord[0][1]),1)[0]
        sizeQ_rd = self.dlg.spinBox_6.value()
        sizeM_rd = self.dlg.spinBox_7.value()
        
        if self.dlg.targetButton.isChecked():
            createTargeted(selectedReliefLayer, vis_lay, px_size, selectedSlopeLayer, target_coord, source_coord, lumen, beamAngle, source_height, sizeQ_rd, sizeM_rd)
        else: 
            if self.dlg.viewButton.isChecked():
                createViewshed(selectedReliefLayer, vis_lay, px_size, source_coord, lumen, beamAngle, source_height, sizeQ_rd, sizeM_rd)
            else:
                createVisibility(selectedReliefLayer,selectedBuildingLayer,vis_lay,source_coord[0],source_height,px_size,sizeQ_rd,sizeM_rd)
        self.dlg.runButton.setEnabled(True)

    def run(self):
        """Run method that performs all the real work"""

        # Create the dialog with elements (after translation) and keep reference
        # Only create GUI ONCE in callback, so that it will only load when the plugin is started
        if self.first_start == True:
            self.first_start = False
            self.dlg = ViewshedQOSDialog()

        # Select a Relief layer 
        #Fetch the current loaded layer
        # layers = QgsProject.instance().layerTreeRoot().children()
        raster_layers = getAllRasterLayers()
        #Clear content
        self.dlg.comboBox.clear()
        #Add all layers
        self.dlg.comboBox.addItems([layer.name() for layer in raster_layers])
        
        # Select a Slope layer 
        #Clear content
        self.dlg.comboBox_5.clear()
        #Add all layers
        self.dlg.comboBox_5.addItems([layer.name() for layer in raster_layers])
        
        # Select a Visibilty layer 
        #Clear content
        self.dlg.comboBox_6.clear()
        #Add all layers
        self.dlg.comboBox_6.addItems([layer.name() for layer in raster_layers])
        
        # Select a Building layer 
        buildingsLayer =  [l.name() for l in QgsProject.instance().mapLayers().values() if l.type() == QgsVectorLayer.VectorLayer]
        #Clear content
        self.dlg.comboBox_7.clear()
        #Add all layers
        self.dlg.comboBox_7.addItems([layer for layer in buildingsLayer])
        
        # Selecting the target
        #Fetch the current loaded layer
        point_layers = getAllPointLayers()
        #Clear content
        self.dlg.comboBox_2.clear()
        #Add all layers
        self.dlg.comboBox_2.addItems([layer.name() for layer in point_layers])
        
        # Selecting an origin
        #Clear content
        self.dlg.comboBox_3.clear()
        #Add all layers 
        self.dlg.comboBox_3.addItems([layer.name() for layer in point_layers])
        
        #Add event handler when the button is clicked
        # self.dlg.toolButton_4.clicked.connect(self.open_hardware_windows)
        
        #Print the result
        # self.dlg.commandLinkButton.clicked.connect(self.print_results)
        
        #Connection between radio button
        # self.dlg.targetButton.toggled.connect(self.deactivateView)
        # self.dlg.viewButton.toggled.connect(self.deactivateTarget)
        
        #Set correct range
        self.dlg.spinBox.setRange(0,100000)
        self.dlg.spinBox_2.setRange(0,1000000000)
        self.dlg.spinBox_3.setRange(0,180)
        self.dlg.spinBox_4.setRange(0,10000)
        self.dlg.spinBox_5.setRange(0,100000)
        self.dlg.spinBox_6.setRange(0,100000)
        self.dlg.spinBox_7.setRange(0,1000000)
        
        #Add event listener
        self.dlg.pushButton_3.clicked.connect(self.addLumenValue)
        self.dlg.pushButton_4.clicked.connect(self.addCandelaValue)
        self.dlg.pushButton_5.clicked.connect(self.addBeamAngleValue)
        self.dlg.runButton.clicked.connect(self.createRaster)
        
        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            
            pass
