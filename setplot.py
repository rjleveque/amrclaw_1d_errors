
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from __future__ import absolute_import
from __future__ import print_function

betal = 5.
al = 1.
ar = 0.
x1 = -10.
x2 = -1.

amr_color = ['g','b','r','c','m','g','b','r','c','m']
amr_marker = ['x','s','o','^','<','.','.','.','.','.']
amr_linestyle = ['-','-','-','-','-','-','-','-','-','-']
amr_plotstyle = [ma+li for ma,li in zip(amr_marker,amr_linestyle)]
print('amr_plotstyle = ',amr_plotstyle)
#amr_plotstyle = ['x-','s-','o-','^-','.-']


#--------------------------
def setplot(plotdata=None):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    plotdata.clearfigures()  # clear any old figures,axes,items data

    def draw_interface_add_legend(current_data):
        from pylab import plot
        from numpy import abs, where, log10, exp, sin, linspace
        #plot([0., 0.], [-1000., 1000.], 'k--')
        try:
            from clawpack.visclaw import legend_tools
            labels = ['Level 1','Level 2', 'Level 3', 'Level 4', 'Level 5',
                      'Level 6','Level 7', 'Level 8', 'Level 9', 'Level 10']
            legend_tools.add_legend(labels, colors=amr_color,
                        markers=amr_marker, linestyles=amr_linestyle,
                        loc='upper left')
        except:
            pass

        # exact solution:
        t = current_data.t
        xx = linspace(-12,12,10000)
        #xpct = xx + t
        #xmct = xx - t
        #p_true = ar*exp(-betar*(xpct-5)**2) * sin(freqr*xpct) + \
        #         al*exp(-betal*(xmct+5)**2) * sin(freql*xmct)
        p_true = p_true_fcn(xx,t)
        plot(xx,p_true,'k')

    def draw_interface_add_legend_innerprod(current_data):
        from pylab import plot
        #plot([0., 0.], [-1000., 1000.], 'k--')
        try:
            from clawpack.visclaw import legend_tools
            labels = ['Level 3','Level 4']
            legend_tools.add_legend(labels, colors=['r','c'],
                                    markers=['o','^'], linestyles=['',''],
                                    loc='upper left')
        except:
            pass

    def add_grid(current_data):
        from pylab import grid
        grid(True)

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Pressure and Velocity', figno=1)
    plotfigure.kwargs = {'figsize': (8,8)}
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'   # top figure
    plotaxes.xlimits = [-12,12]
    plotaxes.ylimits = [-1.1,1.1]
    plotaxes.title = 'Pressure'
    plotaxes.afteraxes = draw_interface_add_legend

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.amr_color = amr_color
    plotitem.amr_plotstyle = amr_plotstyle
    plotitem.amr_data_show = [1,1,1]
    plotitem.amr_kwargs = [{'markersize':5},{'markersize':4},{'markersize':3}]


    # Figure for error

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'   # bottom figure
    plotaxes.xlimits = [-12,12]
    plotaxes.ylimits = [-10,1]
    plotaxes.title = 'log10(Error)'
    plotaxes.afteraxes = add_grid
    
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = log_error
    plotitem.amr_color = amr_color
    plotitem.amr_plotstyle = amr_plotstyle
    plotitem.amr_data_show = [0,0,0,0,1]
    
    
    # Figure for inner product, q[2]
    
    plotfigure = plotdata.new_plotfigure(name='Inner Product', figno=10)
    plotfigure.show = False
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-12,12]
    plotaxes.ylimits = [-15,1]         # use when taking inner product with forward solution
    #plotaxes.ylimits = [-0.01,0.02]    # use when taking inner product with Richardson error
    plotaxes.title = 'log10(Inner Product)'
    plotaxes.afteraxes = draw_interface_add_legend
    
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = plot_innerprod
    plotitem.amr_color = amr_color
    plotitem.amr_plotstyle = amr_plotstyle
    plotitem.amr_data_show = [0,1,1,1,0]
    plotitem.show = True       # show on plot?

    
    # Figure for abs(error) 

    plotfigure = plotdata.new_plotfigure(name='Error', figno=11)
    plotfigure.show = False
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-12,12]
    plotaxes.ylimits = [-15,1]
    plotaxes.title = 'log10(Error)'
    plotaxes.afteraxes = draw_interface_add_legend
    
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = plot_error
    plotitem.amr_color = amr_color
    plotitem.amr_plotstyle = amr_plotstyle
    plotitem.amr_data_show = [1,1,1,1,1]
    plotitem.show = True       # show on plot?

    
    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                                         type='each_gauge')
    plotfigure.clf_each_gauge = True
    plotfigure.kwargs = {'figsize': (10,10)}
                                         
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Pressure'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Velocity'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 1
    plotitem.plotstyle = 'b-'

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

def plot_innerprod(current_data):
    from numpy import abs, where, log10
    abs_ip = abs(current_data.aux[2,:])
    log_ip = where(abs_ip > 1e-15, log10(abs_ip), -15)
    return log_ip
    
def p_true_fcn(x,t):
    from numpy import abs, exp, sin
    xpct = x + t
    xmct = x - t
    p_true = al/(1+exp(-betal*(xmct-x1)) + exp(betal*(xmct-x2))) \
             * sin(0.5*xmct**2)
    return p_true
    
def log_error(current_data):
    from numpy import abs, where, log10, exp, sin
    t = current_data.t
    p = current_data.q[0,:]
    x = current_data.x
    #xpct = x + t
    #xmct = x - t
    #p_true = al/(1+exp(-betal*(x-x1)) + exp(betal*(x-x2))) * sin(0.5d0*xmct**2)
    p_true = p_true_fcn(x,t)
    abs_err = abs(p - p_true)
    log_err = where(abs_err > 1e-15, log10(abs_err), -15)
    return log_err
