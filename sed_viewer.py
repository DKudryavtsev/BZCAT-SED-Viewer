import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
from scipy.optimize import minimize_scalar

PATH = os.path.dirname(os.path.realpath(__file__))
FILE_SEDS_RES = os.path.join(PATH, 'data/seds_residents.csv')
FILE_SEDS_ALL = os.path.join(PATH, 'data/seds_all.csv')


def fit_poly(df, n, seds_res, seds_all):
    """Fitting a polynomial to an object's SED
       n is the row number
    """
    def f(x):
        return -(w[0] + w[1]*x + w[2]*x**2 + w[3]*x**3)
        
    # Constraints
    LOW_BOUND =10
    HIGH_BOUND = 15.5
    LOW_HIGH_BORDER = 12
    
    source = df.iloc[n-1]['BZCAT5 Source name']
    goodbad = df.iloc[n-1]['Correct']
    sed_res = seds_res[seds_res['Source'] == df.iloc[n-1]['Source']]
    sed_all = seds_all[seds_all['Source'] == df.iloc[n-1]['Source']]
    
    # Sometimes coordinates don't match in Dec arcseconds
    if len(sed_res) == 0: 
        sed_res = seds_res[seds_res['Source'].apply(
            lambda x: x[:-1]) == df.iloc[n-1]['Source'][:-1]]
        sed_all = seds_all[seds_all['Source'].apply(
            lambda x: x[:-1]) == df.iloc[n-1]['Source'][:-1]]
        
    # Data both lower and higher LOW_HIGH_BORDER (but lower HIGH_BOUND) 
    # should be present  
    mask_low = sed_res['log_nu'] < LOW_HIGH_BORDER
    mask_high = ((sed_res['log_nu']>LOW_HIGH_BORDER) 
                 & (sed_res['log_nu']<HIGH_BOUND))
    if (len(sed_res[mask_low])<1) or (len(sed_res[mask_high])<1):
        return np.NaN
    
    # Polynomial fit
    sed_pol = sed_res[sed_res['log_nu'] < HIGH_BOUND]
    sed_pol = sed_pol.groupby('log_nu')['log_nufnu'].mean().reset_index()
    nu_max, f_max = np.NaN, np.NaN
    while ((sed_pol.shape[0] > 3) 
           and (sed_pol[(sed_pol['log_nu']>LOW_HIGH_BORDER)
                        & (sed_pol['log_nu']<HIGH_BOUND)].shape[0] > 0)):
        x0 = np.ones(sed_pol.shape[0])
        x1 = np.array(sed_pol['log_nu'])
        x2 = np.array(sed_pol['log_nu']**2)
        x3 = np.array(sed_pol['log_nu']**3)
        X = np.column_stack((x0, x1, x2, x3))
        y = np.array(sed_pol['log_nufnu'])
        w = np.linalg.inv(X.T@X) @ X.T @ y            
            
        # if the polynomial falls at right, then exit
        if (w[1] + 2*w[2]*HIGH_BOUND + 3*w[3]*HIGH_BOUND**2) < 0: 
            res = minimize_scalar(
                f, bounds=(LOW_BOUND, HIGH_BOUND), method='bounded')
            nu_max = res.x  
            f_max = -f(nu_max)
            break 
                               
        # otherwise remove one point at right
        sed_pol = sed_pol[:-1]  
    
    df.loc[n-1, 'Synch_max'] = nu_max
    df.to_csv(file, index=False)   
    
    return source, nu_max, f_max, w, sed_res, sed_all, goodbad


def plot_sed(ax, source, nu_max, f_max, w, sed_res, sed_all, goodbad, n):
    "Plotting a SED with a fitted polynomial"
    ax[0].clear()
    ax[1].clear()
    x = np.linspace(7.5, 15.5, 50)
    y = w[0] + w[1]*x + w[2]*x**2 +w[3]*x**3
    ax[0].errorbar(
        sed_res['log_nu'], sed_res['log_nufnu'], yerr=sed_res['log_err'], fmt='.')
    ax[0].plot(x, y)
    ax[0].vlines(nu_max, ymin=f_max-3, ymax=f_max+1, linestyle='--', color='r')
    ax[1].errorbar(
        sed_all['log_nu'], sed_all['log_nufnu'], yerr=sed_all['log_err'], fmt='.')
    ax[1].plot(x, y)
    ax[1].vlines(nu_max, ymin=f_max-3, ymax=f_max+1, linestyle='--', color='r')    
    plt.suptitle(f'Source: {source}; object No. {n}; nu_max = {nu_max:.3f}; correctness: {goodbad}')
    ax[0].set_title('Resident catalogs')
    ax[1].set_title('All catalogs')
    ax[0].set(
        xlabel='$\\log_{10} \\nu$, [Hz]', 
        ylabel='$\\log_{10} \\nu F_{\\nu}$, [erg cm$^{-2}$ s$^{-1}$]')
    ax[1].set(xlabel='$\\log_{10} \\nu$, [Hz]')
    fig.canvas.draw()
    return


def on_press(event):
    "Hotkeys"
    global n
    #print('press', event.key)
    sys.stdout.flush()
    if event.key == 'enter' or event.key == ' ':
        n += 1
        params = fit_poly(df, n, seds_res, seds_all)
        plot_sed(ax, *params, n)
    elif event.key == 'backspace':
        n -= 1
        params = fit_poly(df, n, seds_res, seds_all)
        plot_sed(ax, *params, n)
    elif event.key == 'b':
        df.loc[n-1, 'Correct'] = False
        df.to_csv(file, index=False)
        n += 1
        params = fit_poly(df, n, seds_res, seds_all)
        plot_sed(ax, *params, n)
    elif event.key == 'g':
        df.loc[n-1, 'Correct'] = True
        df.to_csv(file, index=False)
        n += 1
        params = fit_poly(df, n, seds_res, seds_all)
        plot_sed(ax, *params, n)
    elif event.key == 'c':
        df.loc[n-1, 'Correct'] = np.NaN
        df.to_csv(file, index=False)
        n += 1
        params = fit_poly(df, n, seds_res, seds_all)
        plot_sed(ax, *params, n)
        

# Data readout
if len(sys.argv)!=2 and len(sys.argv)!=3:
    print('\nUsage:\nsed_viewer.py file.csv row_number\n')
    exit()
if len(sys.argv) == 2:
    n = 1
else:
    n = int(sys.argv[2])
file = sys.argv[1]
df = pd.read_csv(file)
df['Synch_max'] = df.get('Synch_max', np.NaN)
df['Correct'] = df.get('Correct', np.NaN)
seds_res = pd.read_csv(FILE_SEDS_RES)
seds_all = pd.read_csv(FILE_SEDS_ALL)

# Recalculation of SED values to logarithmic scale
seds_res['log_nu'] = np.log10(seds_res['Frequency (Hz)'])
seds_res['log_nufnu'] = np.log10(seds_res['Nufnu (erg cm^-2 s^-1)'])
seds_res['log_err'] = (
    np.log10(seds_res['Nufnu (erg cm^-2 s^-1)'] + seds_res['Nufnu_err'])
    - seds_res['log_nufnu'])
seds_all['log_nu'] = np.log10(seds_all['Frequency (Hz)'])
seds_all['log_nufnu'] = np.log10(seds_all['Nufnu (erg cm^-2 s^-1)'])
seds_all['log_err'] = (
    np.log10(seds_all['Nufnu (erg cm^-2 s^-1)'] + seds_all['Nufnu_err'])
    - seds_all['log_nufnu'])

# "Mutual" source names for the datasets
df['Source'] = df.apply(lambda x:
    (x['RA (J2000.0)'][0:2] + x['RA (J2000.0)'][3:5] + x['RA (J2000.0)'][6:8]
     + x['Dec (J2000.0)'][0:3] + x['Dec (J2000.0)'][4:6] + x['Dec (J2000.0)'][7:9]),
    axis=1)

# User information in the terminal
print("""
      Keyboard commands:
      Enter or Space to see the next object,
      Backspace to return one step back,
      'g' to mark as good approximation,
      'b' to mark as bad approximation,
      'c' to clear a mark,
      'f' full screen mode on/off
      'q' to quit 
      """)

# Making a chart
plt.rc('figure', titlesize=18)
plt.rc('axes', titlesize=16)
plt.rc('axes', labelsize=14)
fig, ax = plt.subplots(1, 2, figsize=(16, 5))
fig.canvas.mpl_connect('key_press_event', on_press)
params = fit_poly(df, n, seds_res, seds_all)
plot_sed(ax, *params, n)
plt.show()
