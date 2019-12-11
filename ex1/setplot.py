
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from __future__ import absolute_import
from __future__ import print_function
import gridtools1
import numpy
from clawpack.clawutil.data import ClawData

amrdata = ClawData()
amrdata.read('amr.data',force=True)
probdata = ClawData()
probdata.read('setprob.data',force=True)


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

maxlevels = amrdata.amr_levels_max
amr_color = amr_color[:maxlevels]
amr_marker = amr_marker[:maxlevels]
amr_linestyle = amr_linestyle[:maxlevels]
amr_plotstyle = amr_plotstyle[:maxlevels]

c = ['palegreen','skyblue','salmon','cyan','thistle','yellowgreen',\
     'gold','hotpink','cornflowerblue','lime']
v_levels = numpy.arange(0.5,maxlevels+1,1)
c_levels = c[:maxlevels]

xout = numpy.linspace(-12,12,12000)
tolerance = amrdata.flag_richardson_tol

print('Using %i Levels with tolerance = %.6f' % (maxlevels,tolerance))

xlimits = [-12,1]
ylimits_error = [1e-7, 1.]


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
        
    def color_by_level(current_data):
        from pylab import vstack,contourf,plot,ones,arange,colorbar
        fs = current_data.framesoln
        pout,level = gridtools1.grid_output_1d(fs, 0, xout, return_level=True)
        Xout = vstack((xout,xout))
        Yout = vstack((-1.1*ones(xout.shape), 1.1*ones(xout.shape)))
        L = vstack((level,level))
        contourf(Xout,Yout,L,v_levels,colors=c_levels) 
        cb = colorbar(ticks=range(1,maxlevels+1))
        cb.set_label('AMR Level')
        plot(xout,pout,'k')   
        #import pdb; pdb.set_trace()
        
    def error_color_by_level(current_data):
        from pylab import vstack,contourf,plot,ones,arange,colorbar,\
                          ylim,semilogy
        fs = current_data.framesoln
        t = current_data.t
        pout,level = gridtools1.grid_output_1d(fs, 0, xout, return_level=True)
        err = abs(pout - p_true_fcn(xout, t))
        Xout = vstack((xout,xout))
        Yout = vstack((ylimits_error[0]*ones(xout.shape), 
                       ylimits_error[1]*ones(xout.shape)))
        L = vstack((level,level))
        contourf(Xout,Yout,L,v_levels,colors=c_levels) 
        cb = colorbar(ticks=range(1,maxlevels+1))
        cb.set_label('AMR Level')
        semilogy(xout,err,'k')
        #semilogy(xout,level,'k')
        if tolerance is not None:
            plot(xout,tolerance*ones(xout.shape),'r--')
        
    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Pressure and Velocity', figno=1)
    plotfigure.kwargs = {'figsize': (8,8)}
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'   # top figure
    plotaxes.xlimits = xlimits
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
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [1e-10,1]
    plotaxes.title = 'abs(Error)'
    plotaxes.afteraxes = add_grid
    
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_semilogy')
    plotitem.plot_var = abs_error
    plotitem.amr_color = amr_color
    plotitem.amr_plotstyle = amr_plotstyle
    plotitem.amr_data_show = [1,1,1,1,1]
    
    plotfigure = plotdata.new_plotfigure(name='Pressure and Error', figno=2)
    plotfigure.show = False
    plotfigure.kwargs = {'figsize': (12,8)}
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'   # top figure
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-1.1,1.1]
    plotaxes.title = 'Pressure'
    plotaxes.beforeaxes = color_by_level
    plotaxes.afteraxes = add_grid #draw_interface_add_legend

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = 0
    plotitem.amr_color = amr_color
    plotitem.amr_plotstyle = amr_plotstyle
    plotitem.amr_data_show = [1,1,1]
    plotitem.amr_kwargs = [{'markersize':5},{'markersize':4},{'markersize':3}]


    # Figure for error

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'   # bottom figure
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits_error
    plotaxes.title = 'abs(Error)'
    plotaxes.beforeaxes = error_color_by_level
    plotaxes.afteraxes = add_grid
    
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_semilogy')
    plotitem.show = False
    plotitem.plot_var = abs_error
    plotitem.amr_color = amr_color
    plotitem.amr_plotstyle = amr_plotstyle
    plotitem.amr_data_show = [1,1,1,1,1]
    

    def plot_finest(current_data):
        from pylab import vstack,contourf,plot,ones,arange,colorbar,\
                          xlim,ylim,semilogy,figure,title,clf,subplot,show,draw,\
                          tight_layout,ylabel,grid
        
        fs = current_data.framesoln
        t = current_data.t
        print('+++ plot_finest at t = %.4f' % t)
        pout,level = gridtools1.grid_output_1d(fs, 0, xout, return_level=True)
        err = abs(pout - p_true_fcn(xout, t))
        Xout = vstack((xout,xout))
        L = vstack((level,level))
        figure(3, figsize=(12,8))
        clf()
        
        subplot(311)
        Yout = vstack((-1.1*ones(xout.shape), 1.1*ones(xout.shape)))
        contourf(Xout,Yout,L,v_levels,colors=c_levels) 
        cb = colorbar(ticks=range(1,maxlevels+1))
        cb.set_label('AMR Level')
        plot(xout,pout,'k')
        xlim(xlimits)
        ylim(-1.1,1.1)
        title('Pressure at t = %.4f' % t)
        
        subplot(312)
        Yout = vstack((ylimits_error[0]*ones(xout.shape), 
                       ylimits_error[1]*ones(xout.shape)))
        contourf(Xout,Yout,L,v_levels,colors=c_levels) 
        cb = colorbar(ticks=range(1,maxlevels+1))
        cb.set_label('AMR Level')
        semilogy(xout,err,'k')
        if tolerance is not None:
            plot(xout,tolerance*ones(xout.shape),'r--')
        xlim(xlimits)
        ylim(ylimits_error)
        ylabel('abs(error)')
        grid(True)

        subplot(313)
        Yout = vstack((0*ones(xout.shape), (maxlevels+1)*ones(xout.shape)))
        contourf(Xout,Yout,L,v_levels,colors=c_levels) 
        cb = colorbar(ticks=range(1,maxlevels+1))
        cb.set_label('AMR Level')
        plot(xout,level,'k')
        xlim(xlimits)
        ylim(0,maxlevels+1)
        ylabel('AMR Level')
        tight_layout()
        grid(True)
        draw()
                
    plotfigure = plotdata.new_plotfigure(name='finest', figno=3)
    plotfigure.kwargs = {'figsize': (12,8)}
    plotdata.afterframe = plot_finest
    
    
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
    plotitem.plot_var = 0 #plot_error
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
             * sin(0.25*abs(xmct)**2.5)
    return p_true
    
def abs_error(current_data):
    from numpy import abs, where, log10, exp, sin
    t = current_data.t
    p = current_data.q[0,:]
    x = current_data.x
    p_true = p_true_fcn(x,t)
    abs_err = abs(p - p_true)
    return abs_err
