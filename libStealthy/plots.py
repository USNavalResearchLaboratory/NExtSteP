#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
from scipy.stats import expon, halfnorm, norm, ttest_ind
from matplotlib.widgets import SpanSelector, Button
from abc import ABCMeta, abstractmethod
from math import sqrt

import matplotlib.pyplot as plt
import numpy as np
import warnings
import random


###########################

class BaseCompareIPDsToIPDs(object, metaclass=ABCMeta):
    '''
    Compare IPD sequence with reference IPDs.

    visIPDsVsIPDs method visualizes two different IPD sequences (e.g., a stego one and a referecne cover one) for manual inspection
    '''

    @abstractmethod
    def detectIPDsVsIPDs(self, obsIPDseq, refIPDseq):
        '''
        Returns True if covert communication detected, False otherwise
        '''
        pass

    def visIPDsVsIPDs(self, obsIPDseq, refIPDseq, nbins=100):
        '''
        Visualize distributions of observed and reference IPDs
        '''

        _, (trace, dist,cum) = plt.subplots(3, 1)

        trace.plot([float(x)/(len(obsIPDseq)-1) for x in range(len(obsIPDseq))], obsIPDseq, 'k-', markersize=3.)
        trace.plot([float(x)/(len(refIPDseq)-1) for x in range(len(refIPDseq))], refIPDseq, 'b+', markersize=3.)

        dist.hist(refIPDseq, nbins, normed=True, histtype='stepfilled')
        dist.hist(obsIPDseq, nbins, normed=True, histtype='step')

        cum.hist(refIPDseq, nbins, normed=1, histtype='stepfilled', cumulative=True)
        cum.hist(obsIPDseq, nbins, normed=1, histtype='step', cumulative=True)

        plt.show()


class TrivialVisIPDsVsIPDs(BaseCompareIPDsToIPDs):
    '''
    Class with trivial detector (always returns False) for use in visualization
    '''

    def detectIPDsVsIPDs(self, obsIPDseq, refIPDseq):
        '''
        Returns False (does not attempt to detect covert communication)
        '''
        return False



def selectSpan(y, x = None, title='Press left mouse button and drag select range'):
    '''
    Graphically Select a range from a list or tuple.

    The basic core of this was gorked from this implementation:
    https://matplotlib.org/examples/widgets/span_selector.html

    Args:
        y (list or tuple):
            A list or tuple containing the data we will select from.
        x (list or tuple):
            The x axis of the chart to plot.  If not provided x = range(len(y))
        title (str):
            The title of the chart.

    Raises:
        IndexError:
            Raised if y and x are not the same length.

    Returns:
        (tuple of ints):A tuple containing the starting and stoping
        indexes selcted.

    '''
    # Used to maintain some state info about what is selected so we
    # can update the textbox on the plot.
    class RangeFinder(object):
        boxTxt = 'Idx:\n    Min: %d\n    Max: %d'#\nVal at:\n    Min: %f\n    Max: %f'
        def __init__(self, ax = None, props = None):
            # setup axis
            if ax is None:
                self.ax = plt.axes()
            else:
                self.ax = ax

            # setup box properties
            if props is None:
                self.props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            else:
                self.props = props

            # null out selected range
            self.selectedRange = (0, 0)

            # create text overlay.
            msg = RangeFinder.boxTxt % self.selectedRange
            self.text = self.ax.text(0.05, 0.95, msg, transform=ax.transAxes, verticalalignment='top', bbox=self.props)

        def onSelect(self, xmin, xmax):
            # get indexs from floatin point plot corads
            indmin, indmax = np.searchsorted(x, (xmin, xmax))
            indmax = min(len(x) - 1, indmax)

            # update textbox
            self.selectedRange = (indmin, indmax)
            msg = RangeFinder.boxTxt % self.selectedRange
            self.text.set_text(msg)

    # confirmtion button callback
    def confirm(event):
        plt.close()

    # setup x axis
    if x is None:
        x = range(len(y))
    elif len(x) != len(y):
        raise IndexError('x and y are not of the same length.  x=%d, y=%d' % (len(x), len(y)))

    # setup plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(y, '-')
    ax.set_title(title)

    plt.yscale('log')

    # Setup confirm button
    axConfirm = plt.axes([0.05, 0.9, 0.1, 0.075])
    bnConfirm = Button(axConfirm, 'Confirm')
    bnConfirm.on_clicked(confirm)

    selector = RangeFinder(ax)
    # set useblit True on gtkagg for enhanced performance
    _ = SpanSelector(ax, selector.onSelect, 'horizontal', useblit=True,
                    rectprops=dict(alpha=0.5, facecolor='red'), span_stays=True)
    plt.show()
    # convert from numpy.int64 to python ints
    return tuple(map(int, selector.selectedRange))

def checkOverlap(y, starts, ends, x=None, colors=('r', 'g', 'b', 'y', 'c'), title='Overlap', pdf=False):
    '''
    Create a chart showing the overlap of the provided ranges.

    Args:
        y (list or tuple):
            The data the ranges are over.
        starts (list or tuple):
            The start indexes of the ranges.
        ends (list or tuple):
            The end indexes of the ranges.
        x (list or tuple):
            The x axis of the chart to plot.  If not provided x = range(len(y))
        colors (list or tuple):
            The colors the non-overlaing ranges of the chart will
            cycle through.
        title (str):
            The title the chart will display.
        pdf (bool):
            Should the chart be output as a PDF.  If this is running
            on a headless system this should be set to True.

    Raises:
        IndexError:
            Raised if y and x are not the same length.
    '''
    txtStr = 'Overlap for: %d'
    text = None
    # callback to update message box
    def onpick(event):
        msg = txtStr % int(event.artist.overlap)
        text.set_text(msg)
        plt.draw()
    
    # setup x axis
    if x is None:
        x = range(len(y))
    elif len(x) != len(y):
        raise IndexError('x and y are not of the same length.  x=%d, y=%d' % (len(x), len(y)))

    spans = []

    # start plot setup
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x,y, zorder=-1)
    ax.set_title(title)

    # create overlap rectangles
    zipped = zip(starts, ends) ####
    for i, start, end in zip(range(len(starts)), starts, ends):
        overlap = 0
        if i < len(ends)-1:
            nextStart = starts[i+1]
            # check to see if we overlap with the start of the next
            # range provided.
            if end > nextStart:
                # calc the bounds of the hatched rectangle
                overlap = end - nextStart
                print("-----------Overlap ", i, " is ", overlap) #### 
                newEnd = end - overlap
                print("--------- New End ", i, " is ", newEnd)  #### 

                # create hatched overlap rectangle
                rectprops = dict(fill=False, hatch='/', picker=True, zorder=i+1)
                poly = plt.axvspan(newEnd, end, **rectprops)
                # piggyback on the polygon object so that we can pass
                # the overlap value to our text update callback
                poly.overlap = overlap

                spans.append(poly)
                end = newEnd

        # Draw non-overlaping rectangle.
        rectprops = dict(alpha=0.25, color=colors[i%len(colors)], zorder=i,picker=True)
        poly = plt.axvspan(start, end, **rectprops)
        # piggyback on the polygon object so that we can pass
        # the overlap value to our text update callback
        poly.overlap = overlap

        spans.append(poly)

    plt.yscale('log')
    # textbox properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # setup the default message for the textbox.
    largestOverlap = max(spans, key=lambda x: x.overlap).overlap
    msg = txtStr % largestOverlap
    print("---------------Largest Overlap: ", msg) #### 
    text = ax.text(0.05, 0.95, msg, transform=ax.transAxes, verticalalignment='top', bbox=props)

    if pdf:
        plt.savefig(pdf)
    else:
        fig.canvas.mpl_connect('pick_event', onpick)
        plt.show()
