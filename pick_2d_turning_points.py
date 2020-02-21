import numpy as np
from pathlib import Path
from astropy.table import Table, vstack, unique
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.time import Time
from datetime import datetime
import j0600

band = 'rp'
#  B
#   I
#   R
#   V
#  gp
#  ip
#  rp
fig, ax = plt.subplots()

print('press q to quit')

inf = 'turning_points.txt'
outf = 'turning_points_{}.txt'.format(band)
print('Input file is {}'.format(inf))
print(outf)
my_file = Path(inf)
if my_file.is_file():
    print('{} exists! Reading in...'.format(my_file))
    t = Table.read(inf,format='ascii.no_header')
    xlines=t['col1']
    ylines=t['col2']
    print(t)
else:
    xlines=np.array([58880.])
    ylines=np.array([14.])
    
# read in all the data!
(t, dm) = j0600.j0600()

# reject ones with large errors
t = t[t['Uncertainty']<0.06]

tg = t[np.where(t['Band']==band)]


print(xlines)

class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration
    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args: # changed for v1.5, see below
            self.fall = True
            return True
        else:
            return False


def onpress(event):
    global xlines, ylines
    global ax
    #print('you pressed {} at position X={} Y={}'.format(event.key, event.xdata, event.ydata))

    for case in switch(event.key):

        ###print('key pressed!')
        if (event.xdata==None):
            break
        if (event.ydata==None):
            break
        
        idx = np.abs(xlines - event.xdata).argmin()
        if case('a'):
            # add line at that location
            xlines = np.append(xlines, event.xdata)
            ylines = np.append(ylines, event.ydata)

            break
        
        if case('m'):
            # find nearest value and update that
            xlines[idx] = event.xdata
            ylines[idx] = event.ydata
            break
        
        if case('y'):
            # find nearest value and update that
            ylines[idx] = event.ydata
            break

        if case('d'):
             # find nearest value and delete that
            xlines = np.delete(xlines, idx)
            ylines = np.delete(ylines, idx)

        if case('w'):
            # write out the table
            
            t=Table([xlines, ylines])
            t.write(outf, format='ascii.no_header', overwrite=True)
            print('Written {} out!'.format(outf))
            print(t)

        if case('q'):
            print('q is for quitters')
            break

        if case():
            print('not a recognised keypress')
            break

    # set them up in an ordered list
    s_ind = np.argsort(xlines)
    xlines = xlines[s_ind]
    ylines = ylines[s_ind]
    
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()

    ax.clear()
    print(xlines)
    
    ax.vlines(xlines,0,20)
    ax.errorbar(tg['MJD'],tg['Magcorr'], yerr=tg['Uncertainty'],fmt='b.')
    ax.plot(xlines,ylines,'r')


    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    plt.draw()


ax.errorbar(tg['MJD'],tg['Magcorr'], yerr=tg['Uncertainty'],fmt='b.')

ax.vlines(xlines,0,20)
ax.plot(xlines,ylines,'r')

cid = fig.canvas.mpl_connect('key_press_event', onpress)

plt.show()

