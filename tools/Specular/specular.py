import wx
import time
import sys
import numpy
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import wx.lib.inspection
#from wx.lib.mixins.listctrl import ColumnSorterMixin
from wx.lib.mixins.listctrl import ListCtrlAutoWidthMixin
import challa.workspace as ws
from astropy.nddata import NDData
from astropy.table import Table
from matplotlib import pyplot as plt  

def load_test():
   ws.import_file("fits/calibrated.ms.image.spectrum.J113740.6-010454.spw0.image.fits")
   ws.import_file("fits/calibrated.ms.line.spw0.source15.image.fits")
   ws.import_file("fits/NGC6240_continuum.fits")
   ws.import_file("fits/logfile_alma_hatlas_cycle1_inc-z_beye.fits")

class AutoWidthListCtrl(wx.ListCtrl, ListCtrlAutoWidthMixin):
    def __init__(self, parent):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT|wx.LC_ALIGN_LEFT|wx.LC_HRULES|wx.LC_VRULES)
        ListCtrlAutoWidthMixin.__init__(self)


class Specular(wx.Frame):
   def __init__(self,parent,title):
      wx.Frame.__init__(self,parent,title=title,size=(1024,768))
      self.createPanels()
      self.createMenuBar()
      self.createToolBar()

   def createToolBar(self):
      toolbar = self.CreateToolBar()
      qtool = toolbar.AddLabelTool(wx.ID_ANY, 'Quit', wx.Bitmap('quit.png'))
      toolbar.Realize()
      self.Bind(wx.EVT_TOOL, self.OnQuit, qtool)

   def createMenuBar(self):
      menubar = wx.MenuBar()
      fileMenu = wx.Menu()
      fitem = fileMenu.Append(wx.ID_EXIT, 'Quit', 'Quit application')
      menubar.Append(fileMenu, '&File')
      self.SetMenuBar(menubar)
      self.Bind(wx.EVT_MENU, self.OnQuit, fitem)

   def createFigures(self):
      self.figure1 = matplotlib.figure.Figure()
      self.figure2 = matplotlib.figure.Figure()
      self.figure3 = matplotlib.figure.Figure()
      self.figure4 = matplotlib.figure.Figure()
      self.mainPlot = self.figure1.add_subplot(111)
      self.auxPlot  = self.figure2.add_subplot(111)
      self.analysisPlot  = self.figure3.add_subplot(111)
      self.resultPlot  = self.figure4.add_subplot(111)
   
   def updateWorkspace(self):
      items = ws.elements().items()
      self.myRowDict=dict()
      for key, data in items:
         print key
         if isinstance(data,Table):
            index = self.wslist.InsertStringItem(sys.maxint, key)
            self.wslist.SetStringItem(index, 1, "Table")
            self.wslist.SetStringItem(index, 2, "---")
            self.myRowDict[index] = key
         elif  isinstance(data,NDData):
            (dim,shape,otype)=ws.real_dims(data)
            index = self.wslist.InsertStringItem(sys.maxint, key)
            self.wslist.SetStringItem(index, 1, otype)
            self.wslist.SetStringItem(index, 2, str(shape))
            self.myRowDict[index] = key

   def createPanels(self):
      self.createFigures()
      
      #Create main panels
      self.primary_panel = wx.Panel(self)
      main_box = wx.BoxSizer(wx.HORIZONTAL)
      self.mainPanel=wx.Panel(self)
      inner_box = wx.BoxSizer(wx.VERTICAL)
     
      stc=wx.StaticBox(self.mainPanel,label="Main View")
      self.mainImgPanel=FigureCanvas(self.mainPanel,-1,self.figure1)
      self.mainImgPanel.SetMinSize((10,10))
      self.sld = wx.Slider(self.mainPanel, value=50, minValue=0, maxValue=100,style=wx.SL_HORIZONTAL|wx.SL_LABELS)
      self.sld.Enable(False)
      self.sld.Bind(wx.EVT_SCROLL, self.OnSliderScroll)
      inner_box.Add(stc)
      inner_box.Add(self.mainImgPanel,flag=wx.EXPAND|wx.TOP,proportion=1)
      inner_box.Add(self.sld,flag=wx.EXPAND|wx.BOTTOM)
      self.mainPanel.SetSizer(inner_box)

      self.auxPanel=FigureCanvas(self.primary_panel,-1,self.figure2)
      self.auxPanel.SetMinSize((10,10))

      #Create list of objects
      self.workspace = wx.Panel(self)
      inner_box = wx.BoxSizer(wx.VERTICAL)
      
      self.wslist = AutoWidthListCtrl(self.primary_panel)
      self.wslist.InsertColumn(0, 'cubename')#, width=500)
      self.wslist.InsertColumn(1, 'type')#, width=40)
      self.wslist.InsertColumn(2, 'dims')#, wx.LIST_FORMAT_RIGHT, 30)
      self.updateWorkspace()
      
      #Dynamics of the panel
      self.wslist.Bind(wx.EVT_LIST_ITEM_SELECTED, self.onItemSelected)

      #Panel Layout
      main_box.Add(self.wslist,flag=wx.EXPAND|wx.LEFT,proportion=1,border=10)
      main_box.Add(self.mainPanel,flag=wx.EXPAND,proportion=1,border=10)
      main_box.Add(self.auxPanel,flag=wx.EXPAND|wx.RIGHT,proportion=1,border=10)
      self.primary_panel.SetSizer(main_box)
      
      #Create analysis panels
      self.secondary_panel = wx.Panel(self, 1)
      second_box = wx.BoxSizer(wx.HORIZONTAL)
      self.analysisPanel=FigureCanvas(self.secondary_panel,-1,self.figure3)
      self.analysisPanel.SetMinSize((10,10))
      self.resultPanel=FigureCanvas(self.secondary_panel,-1,self.figure3)
      self.resultPanel.SetMinSize((10,10))
      self.cpListBox = wx.ListBox(self.secondary_panel, -1)
      
      #Panel Layout
      second_box.Add(self.cpListBox,flag=wx.EXPAND|wx.LEFT,proportion=1,border=10)
      second_box.Add(self.analysisPanel,flag=wx.EXPAND,proportion=1,border=10)
      second_box.Add(self.resultPanel,flag=wx.EXPAND|wx.RIGHT,proportion=1,border=10)
      self.secondary_panel.SetSizer(second_box)

      #General Layout
      frame_box = wx.BoxSizer(wx.VERTICAL)
      frame_box.Add(self.primary_panel,flag= wx.EXPAND|wx.TOP,proportion=1,border=10)
      frame_box.Add(self.secondary_panel,flag= wx.EXPAND|wx.BOTTOM,proportion=1,border=10)
      self.SetSizer(frame_box)
      
      #Status Bar
      self.statusbar = self.CreateStatusBar()
      self.statusbar.SetStatusText("Ready.")

   def PlotMain(self,ndd):
      (dim,shape,otype)=ws.real_dims(ndd)
      self.sld.SetValue(50)
      self.sld.SetMax(100)
      if dim == 1:
         self.mainPlot.plot(ndd)
         self.sld.Enable(False)
      elif dim == 2:
         img=ndd
         if ndd.ndim==4:
            img=ndd[0][0]
         if ndd.ndim==3:
            img=ndd[0]
         self.mainPlot.imshow(img)
         self.sld.Enable(False)
      elif dim == 3:
         if ndd.ndim==4:
            self.target=ndd[0]
         self.sld.SetMax(shape[0]-1)
         self.sldIndex=int(shape[0]/2)
         self.sld.SetValue(self.sldIndex)
         self.sld.Enable(True)
         self.mainPlot.imshow(self.target[self.sldIndex])
 
   def OnSliderScroll(self, e):
      obj = e.GetEventObject()
      self.sldIndex = obj.GetValue()
      self.mainPlot.clear()
      self.mainPlot.imshow(self.target[self.sldIndex])
      self.mainImgPanel.draw()
      self.mainImgPanel.Refresh()


   def onItemSelected(self,event):
      currentItem = event.m_itemIndex
      key=self.myRowDict[currentItem]
      print key
      items=ws.elements()
      self.mainPlot.clear()
      self.PlotMain(items[key])
      self.mainImgPanel.draw()
      self.mainImgPanel.Refresh()

 
   def OnQuit(self, e):
      self.Close()

   def Load(self,event):
     self.statusbar.SetStatusText("Loading...")
     time.sleep(3)
     self.statusbar.SetStatusText("Ready.")

load_test()
app = wx.App(redirect=False)
frame=Specular(None,"SpeCuLar")
frame.Show()
frame.Centre()
wx.lib.inspection.InspectionTool().Show()
app.MainLoop()
