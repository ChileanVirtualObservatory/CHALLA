import wx
import time
import numpy
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

class CubeImagePanel(wx.Panel):
   def __init__(self,parent,size):
      wx.Panel.__init__(self,parent,-1,size)
      self.figure = matplotlib.figure.Figure()
      self.axes = self.figure.add_subplot(111)
      t= numpy.arange(0.0,10,1.0)
      s = [0,1,0,1,0,2,1,2,1,0]
      self.y_max=10
      self.axes.plot(t,s)
      self.canvas = FigureCanvas(self,-1,self.figure)

class Specular(wx.Frame):
   def __init__(self,parent,title):
      wx.Frame.__init__(self,parent,title=title,size=(1024,768))
      #text = wx.StaticText(self,label="Hello, World")
      self.sp = wx.SplitterWindow(self)
      self.appsplit = wx.SplitterWindow(self.sp)
      self.toolpanel = wx.Panel(self.sp, style=wx.SUNKEN_BORDER)
      self.sp.SplitVertically(self.appsplit,self.toolpanel,800)
      self.cubesplit= wx.SplitterWindow(self.appsplit)
      self.clumpsplit= wx.SplitterWindow(self.appsplit)
      self.appsplit.SplitHorizontally(self.cubesplit,self.clumpsplit,384)
      self.cubeimagePanel = CubeImagePanel(self.cubesplit,size=(300,300))
      self.workspacePanel = wx.Panel(self.cubesplit, style=wx.SUNKEN_BORDER)
      self.cubesplit.SplitVertically(self.workspacePanel,self.cubeimagePanel,200)
      self.statusbar = self.CreateStatusBar()
      self.statusbar.SetStatusText("Ready.")

      self.loadButton = wx.Button(self.toolpanel, -1, "Load",size= (100,40), pos=(10,10))
      self.loadButton.Bind(wx.EVT_BUTTON,self.Load)
   
   def Load(self,event):
     self.statusbar.SetStatusText("Loading.")
     time.sleep(3)
     self.statusbar.SetStatusText("Ready.")

app = wx.App(redirect=False)
frame=Specular(None,"SpeCuLar")
frame.Show()
app.MainLoop()
