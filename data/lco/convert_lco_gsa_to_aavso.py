import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table, vstack, unique
from pathlib import Path

# read in aavso table as a template
d = ascii.read('00demo_aavso.txt')

outf = 'obs_LCOGT_GSA.ecsv'

# get a list of all the AAVSO files in the directory

t = ascii.read('00demo_aavso.txt')
t['lco_filename']=['verylongfilenameusedtomakestaaasduffgoon']
t['MJD'] = 12345.00
t['LAT'] = 0.00
t['LON'] = 0.00
t.replace_column('Band',['3lf'])
t.replace_column('Star Name',['ASASSN-V J060000.76-310027.83'])
t.replace_column('Comp Star 1',['UC4 296-009008'])
t.remove_column('Charts')
t.remove_row(0)

# process all the files in the list and add them to the table
# grep all of gary saccos text files
for f in Path('./').rglob('*KB*.txt'):
    
    print(f.name)

    band = f.name[0]
    obs = f.name[2:5]
    print('band is {} and obs is {}'.format(band, obs))

    # convert input table to columns in template table
    tin = ascii.read(f)

    tin.remove_column('col1')
    tin.remove_column('col6')
    tin.remove_column('col5')
    tin.rename_column('col2','lco_filename')
    tin.rename_column('col4','relfluxT1n')
    tin.rename_column('col7','relfluxT1nerr')
    
    tin['MJD'] = tin['col3'] - 0.5
    tin.remove_column('col3')
    tin['Observer Affiliation'] = obs
    tin['Band'] = band
    tin['Observer Code'] = 'GSA'

    tin['Star Name'] = 'ASASSN-V J060000.76-310027.83'

    tin['Magnitude'] = 13 - 2.5 * np.log10(tin['relfluxT1n'])

    ta = 13 - 2.5 * np.log10(tin['relfluxT1n']+tin['relfluxT1nerr'])
    tb = 13 - 2.5 * np.log10(tin['relfluxT1n']-tin['relfluxT1nerr'])
    tin['Uncertainty'] = np.abs(ta-tb)

    tin.remove_column('relfluxT1n')
    tin.remove_column('relfluxT1nerr')

    t = vstack([t,tin])


t = unique(t, 'MJD')

t.sort('MJD')


# write out the table
t.write(outf, format='ascii.ecsv', overwrite=True)
print(t)
print('wrote out table to {}'.format(outf))


