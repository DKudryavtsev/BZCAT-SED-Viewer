import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
import sys
import os
import warnings
from scipy.optimize import minimize_scalar

PATH = os.path.dirname(os.path.realpath(__file__))
FILE_SEDS_RES = os.path.join(PATH, 'data/seds_residents.csv')
FILE_SEDS_ALL = os.path.join(PATH, 'data/seds_all.csv')
LOW_OPT_BOUND = 10    # Left frequency bound (for the optimizer only)
LOW_HIGH_BORDER = 12  # Fluxes lower and higher this freq. must be present
LEFT_PLOT_BOUND = 7.5 # Left bound for polynomial plotting
LOGNUFNU_MAX = -5     # To remove flux outliers


def check_points(sed, degree, high_bound):
    "Checking if the number of points is sufficcient"
    mask1 = sed.shape[0] > degree
    mask2 = (sed[(sed['log_nu']>LOW_HIGH_BORDER) & (sed['log_nu']<high_bound)]
             .shape[0] > 0)
    mask3 = sed[(sed['log_nu'] < LOW_HIGH_BORDER)].shape[0] > 0
    return mask1 and mask2 and mask3


def fit3deg(sed_pol, high_bound):
    "Fitting a 3rd degree polynomial"
    def f(x):
        return -(w[0] + w[1]*x + w[2]*x**2 + w[3]*x**3)
    
    nu_max, f_max, w = np.NaN, np.NaN, [np.NaN, np.NaN, np.NaN, np.NaN]
    while check_points(sed_pol, 3, high_bound):
        x0 = np.ones(sed_pol.shape[0])
        x1 = np.array(sed_pol['log_nu'])
        x2 = np.array(sed_pol['log_nu']**2)
        x3 = np.array(sed_pol['log_nu']**3)
        X = np.column_stack((x0, x1, x2, x3))
        y = np.array(sed_pol['log_nufnu'])
        w = np.linalg.inv(X.T@X) @ X.T @ y            
            
        # if the polynomial falls at right, then exit
        if (w[1] + 2*w[2]*high_bound + 3*w[3]*high_bound**2) < 0: 
            res = minimize_scalar(
                f, bounds=(LOW_OPT_BOUND, high_bound), method='bounded')
            nu_max = res.x  
            f_max = -f(nu_max)
            break 
                               
        # otherwise remove one point at right
        sed_pol = sed_pol[:-1]
        
    return nu_max, f_max, w


def fit2deg(sed_pol, high_bound):
    "Fitting a 2nd degree polynomial"
    def f(x):
        return -(w[0] + w[1]*x + w[2]*x**2)
    
    nu_max, f_max, w = np.NaN, np.NaN, [np.NaN, np.NaN, np.NaN]
    while check_points(sed_pol, 2, high_bound):
        x0 = np.ones(sed_pol.shape[0])
        x1 = np.array(sed_pol['log_nu'])
        x2 = np.array(sed_pol['log_nu']**2)
        X = np.column_stack((x0, x1, x2))
        y = np.array(sed_pol['log_nufnu'])
        w = np.linalg.inv(X.T@X) @ X.T @ y            
            
        # if the polynomial falls at right, then exit
        if (w[1] + 2*w[2]*high_bound) < 0: 
            nu_max = -w[1] / (2*w[2])  
            f_max = -f(nu_max)
            break 
                               
        # otherwise remove one point at right
        sed_pol = sed_pol[:-1]
        
    return nu_max, f_max, w


def fit_poly(df, n, seds_res, seds_all, high_bound, polydeg):
    """Fitting a polynomial to an object's SED
       n is the row number in the table (df)
       
       The function prepares the data and calls fit2deg or fit3deg 
    """
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
       
    # Polynomial fit
    sed_pol = sed_res[sed_res['log_nu'] < high_bound]
    sed_pol = sed_pol.groupby('log_nu')['log_nufnu'].mean().reset_index()
    if polydeg == 3:  
        nu_max, f_max, w = fit3deg(sed_pol, high_bound)
    else:             
        nu_max, f_max, w = fit2deg(sed_pol, high_bound)
    
    df.loc[n-1, 'Synch_max'] = nu_max
    df.loc[n-1, 'High_bound'] = high_bound
    df.loc[n-1, 'Poly_degree'] = polydeg
    df.to_csv(file, index=False)   
    
    return source, nu_max, f_max, w, sed_res, sed_all, goodbad


def plot_sed(ax, source, nu_max, f_max, w, sed_res, sed_all, 
             goodbad, n, high_bound, polydeg):
    "Plotting a SED with a fitted polynomial"
    ax[0].clear()
    ax[1].clear()
    ax[0].errorbar(
        sed_res['log_nu'], sed_res['log_nufnu'], yerr=sed_res['log_err'], fmt='.')
    ax[1].errorbar(
        sed_all['log_nu'], sed_all['log_nufnu'], yerr=sed_all['log_err'], fmt='.')
        
    if np.isnan(w[0]):
        warnings.warn('Insufficient data')
    else:
        x = np.linspace(LEFT_PLOT_BOUND, high_bound, 50)
        if polydeg == 3:
            y = w[0] + w[1]*x + w[2]*x**2 +w[3]*x**3
        else:
            y = w[0] + w[1]*x + w[2]*x**2
        ax[0].plot(x, y)
        ax[0].vlines(nu_max, ymin=f_max-3, ymax=f_max+1, linestyle='--', color='r')
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


def add_n(val):
    "Add/subtract the row number and check if it's within the range"
    global n
    n += val
    if n<1 or n>3561:
        n -= val
        warnings.warn('Object No. is out of range, must be [1:3561]')


def on_press(event):
    "Hotkeys"
    global n
    #print('press', event.key)
    sys.stdout.flush()
    if event.key == 'enter' or event.key == ' ':
        add_n(1)
        update(True)
    elif event.key == 'backspace':
        add_n(-1)
        update(True)
    elif event.key == 'b':
        df.loc[n-1, 'Correct'] = False
        df.to_csv(file, index=False)
        update(True)
    elif event.key == 'g':
        df.loc[n-1, 'Correct'] = True
        df.to_csv(file, index=False)
        update(True)
    elif event.key == 'c':
        df.loc[n-1, 'Correct'] = np.NaN
        df.to_csv(file, index=False)
        update(True)


def get_polydeg(label):
    "Get the polynomial degree from the button labels"
    poly_dict = {'2nd degree': 2, '3rd degree': 3}
    return poly_dict[label]


def set_true_values():
    "Setting saved values from the dadaframe"
    hbslider_val = df.iloc[n-1]['High_bound']
    hb_slider.eventson = False
    hb_slider.set_val(hbslider_val)
    hb_slider.eventson = True
      
    position_dict = {2:0, 3:1}
    polydeg = int(df.iloc[n-1]['Poly_degree'])
    button.eventson = False
    button.set_active(position_dict[polydeg])
    button.eventson = True
    
    return hbslider_val, polydeg


def update(val):
    "Response if the slider or button have been changed by the user"
    if not np.isnan(df.iloc[n-1]['Correct']):
        hbslider_val, polydeg = set_true_values()        
    else:
        hbslider_val = hb_slider.val
        polydeg = get_polydeg(button.value_selected)
    params = fit_poly(df, n, seds_res, seds_all, hbslider_val, polydeg)
    plot_sed(ax, *params, n, hbslider_val, polydeg)


def read_seds(file):
    "Reading SEDs from a file"
    seds = pd.read_csv(file)
 
    # Recalculation of SED values to logarithmic scale
    seds['log_nu'] = np.log10(seds['Frequency (Hz)'])
    seds['log_nufnu'] = np.log10(seds['Nufnu (erg cm^-2 s^-1)'])
    seds['log_err'] = (
        np.log10(seds['Nufnu (erg cm^-2 s^-1)'] + seds['Nufnu_err'])
        - seds['log_nufnu'])
    
    # Removing flux ouliers
    seds = seds[seds['log_nufnu'] < LOGNUFNU_MAX]
    
    return seds

   
# MAIN PROGRAM
# Data readout
if len(sys.argv)!=2 and len(sys.argv)!=3:
    print('\nUsage:\nsed_viewer.py file.csv row_number\n')
    exit()
if len(sys.argv) == 2:
    n = 1
else:
    n = int(sys.argv[2])
    if n<1 or n>3561:
        warnings.warn('Object No. is out of range, must be [1:3561]')
        exit(0)
file = sys.argv[1]

df = pd.read_csv(file)
df['Synch_max'] = df.get('Synch_max', np.NaN)
df['Correct'] = df.get('Correct', np.NaN)
df['High_bound'] = df.get('High_bound', np.NaN)
df['Poly_degree'] = df.get('Poly_degree', np.NaN)

seds_res = read_seds(FILE_SEDS_RES)
seds_all = read_seds(FILE_SEDS_ALL)

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

# adjust the main plot to make room for the slider and button
fig.subplots_adjust(bottom=0.3)

# Make a horizontal slider to control the right frequency bound
axfreq = fig.add_axes([0.2, 0.1, 0.25, 0.03])
hb_slider = Slider(
    ax=axfreq,
    label='Right frequency bound',
    valmin=12.5,
    valmax=20,
    valinit=15.5
)
hb_slider.on_changed(update)

# Radio button: 2nd or 3rd degree polynomial
#axcolor = 'lightgoldenrodyellow'
rax = fig.add_axes([0.5, 0.07, 0.08, 0.1])#, facecolor=axcolor)
button = RadioButtons(rax, ('2nd degree', '3rd degree'), active=1)
button.on_clicked(update)

# Connect the figure to keyboard hotkeys
fig.canvas.mpl_connect('key_press_event', on_press)
update(True)

plt.show()
