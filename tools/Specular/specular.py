import wx
import time
import sys
import numpy
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import wx.lib.inspection
from wx.lib.mixins.listctrl import ColumnSorterMixin
import challa.workspace as ws
from astropy.nddata import NDData
from astropy.table import Table


def load_test():
   ws.import_file("fits/calibrated.ms.image.spectrum.J113740.6-010454.spw0.image.fits")
   ws.import_file("fits/calibrated.ms.line.spw0.source15.image.fits")
   ws.import_file("fits/NGC6240_continuum.fits")
   ws.import_file("fits/logfile_alma_hatlas_cycle1_inc-z_beye.fits")

actresses = {
1 : ('jessica alba', 'pomona', '1981'), 
2 : ('sigourney weaver', 'new york', '1949'),
3 : ('angelina jolie', 'los angeles', '1975'), 
4 : ('natalie portman', 'jerusalem', '1981'),
5 : ('rachel weiss', 'london', '1971'), 
6 : ('scarlett johansson', 'new york', '1984') 
}

class SortedListCtrl(wx.ListCtrl, ColumnSorterMixin):
    def __init__(self, parent):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        ColumnSorterMixin.__init__(self, 30)
        self.itemDataMap = actresses

    def GetListCtrl(self):
        return self


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

   def createPanels(self):
      self.main_panel = wx.Panel(self, 1)

      main_box = wx.BoxSizer(wx.HORIZONTAL)
      self.createFigures()
      self.cubeImagePanel=FigureCanvas(self.main_panel,-1,self.figureA)
      self.cubeImagePanel.SetMinSize((10,10))
      self.secondImagePanel=FigureCanvas(self.main_panel,-1,self.figureA)
      self.secondImagePanel.SetMinSize((10,10))
      #self.wslist = wx.ListCtrl(main_panel, -1, style=wx.LC_REPORT)
      self.wslist = SortedListCtrl(self.main_panel)
      self.wslist.InsertColumn(0, 'cubename', width=100)
      self.wslist.InsertColumn(1, 'type', width=100)
      self.wslist.InsertColumn(2, 'dims', wx.LIST_FORMAT_RIGHT, 90)
      items = ws.elements().items()
      count=0
      for key, data in items:
         if isinstance(data,Table):
            index = self.wslist.InsertStringItem(sys.maxint, key)
            self.wslist.SetStringItem(index, 1, "Table")
            self.wslist.SetStringItem(index, 2, "---")
            self.wslist.SetItemData(index, count)
         elif  isinstance(data,NDData):
            index = self.wslist.InsertStringItem(sys.maxint, key)
            if data.ndim==1:
               self.wslist.SetStringItem(index, 1, "Spectra")
            elif data.ndim==2:
               self.wslist.SetStringItem(index, 1, "Image")
            elif data.ndim==3:
               self.wslist.SetStringItem(index, 1, "Cube")
            else:
               continue
            self.wslist.SetStringItem(index, 2, str(data.shape))
            self.wslist.SetItemData(index, count)
         else:
            continue
         count=count+1

      main_box.Add(self.wslist,flag=wx.EXPAND|wx.LEFT,proportion=1,border=10)
      main_box.Add(self.cubeImagePanel,flag=wx.EXPAND,proportion=1,border=10)
      main_box.Add(self.secondImagePanel,flag=wx.EXPAND|wx.RIGHT,proportion=1,border=10)
      self.main_panel.SetSizer(main_box)

      
      self.second_panel = wx.Panel(self, 1)
      main_box = wx.BoxSizer(wx.HORIZONTAL)
      #self.createFigures()
      self.clumpsImagePanel=FigureCanvas(self.second_panel,-1,self.figureA)
      self.clumpsImagePanel.SetMinSize((10,10))
      self.resultImagePanel=FigureCanvas(self.second_panel,-1,self.figureA)
      self.resultImagePanel.SetMinSize((10,10))
      self.cpListBox = wx.ListBox(self.second_panel, -1)
      main_box.Add(self.cpListBox,flag=wx.EXPAND|wx.LEFT,proportion=1,border=10)
      main_box.Add(self.clumpsImagePanel,flag=wx.EXPAND,proportion=1,border=10)
      main_box.Add(self.resultImagePanel,flag=wx.EXPAND|wx.RIGHT,proportion=1,border=10)
      self.second_panel.SetSizer(main_box)


      frame_box = wx.BoxSizer(wx.VERTICAL)
      frame_box.Add(self.main_panel,flag= wx.EXPAND|wx.TOP,proportion=1,border=10)
      frame_box.Add(self.second_panel,flag= wx.EXPAND|wx.BOTTOM,proportion=1,border=10)

      self.SetSizer(frame_box)

      self.statusbar = self.CreateStatusBar()
      self.statusbar.SetStatusText("Ready.")

   def createFigures(self):
      self.figureA = matplotlib.figure.Figure()
      self.axesA = self.figureA.add_subplot(111)
      #t= numpy.arange(0.0,10,1.0)
      #s = [0,1,0,1,0,2,1,2,1,0]
      #self.y_max=10
      #self.axesA.plot(t,s)
 
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
