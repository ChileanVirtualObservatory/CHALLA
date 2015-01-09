import wx
import matplotlib.figure as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import wx.lib.inspection

class Test(wx.Frame):

    def __init__(self):
        super(Test, self).__init__(parent=None, id=-1)
        self.panel = wx.Panel(self, 1)
        self.panel.SetBackgroundColour('RED')

        self.figure = plt.Figure()
        self.axes1 = self.figure.add_subplot(111)
        self.figurepanel = FigureCanvas(self.panel, -1, self.figure)

        main_box = wx.BoxSizer(wx.HORIZONTAL)
        main_box.Add(self.figurepanel, flag=wx.EXPAND, proportion=1)
        self.panel.SetSizer(main_box)

        frame_box = wx.BoxSizer(wx.VERTICAL)
        frame_box.Add(self.panel, flag=wx.EXPAND, proportion=1)

        self.SetSizer(frame_box)
        self.Show()
        self.Layout()

def main():
    app = wx.App()
    Test()
    wx.lib.inspection.InspectionTool().Show()
    app.MainLoop()

if __name__ == '__main__':
    main()
