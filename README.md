# BZCAT SED Viewer

![](./data/Example.png)

## Description

The BZCAT SED Viewer script displays the spectral energy distributions (SEDs) for the blazars from the [Roma-BZCAT catalog](https://heasarc.gsfc.nasa.gov/W3Browse/all/romabzcat.html) (5th edition, Massaro et al., 2015, Ap&SS, 357, 75). It also calculates the frequency of the synchrotron component maximum with a third-degree polinomial and allows a user to interactively mark it as 'good' or 'bad'. The results are saved in the original file. One can use the total list of Roma-BZCAT objects (data/BZCAT.csv provided here) or make their own subsample (the file must contain the 'BZCAT5 Source name', 'RA (J2000.0)', and 'Dec (J2000.0)' columns).

## Data

* [List of blazars](./data/BZCAT.csv) from the [Roma-BZCAT catalog](https://heasarc.gsfc.nasa.gov/W3Browse/all/romabzcat.html)
* Spectral energy distributions from [SED Builder](https://tools.ssdc.asi.it/SED/)

## Dependencies

Python3 with the numpy, pandas, matplotlib, and scipy libraries.

Installation of the libraries: pip install -r requirements.txt

## Running the script

python3 sed_viewer.py data/BZCAT.csv 1

Here
* data/BZCAT.csv - the provided file with BZCAT objects and their coordinates (or a file with your subsample)
* 1 - the number of the object your want to begin with (the first object in the list corresponds to 1)

The script makes two charts - for the data from the SED Builder resident catalogs and for all catalogs. The 3rd degree polynomial fits the resident catalogs only.

Hotkeys in the graphics window:
* Enter or Space - see the next object,
* Backspace - return one step back,
* 'g' - mark as a good fit,
* 'b' - mark as a bad fit,
* 'c' - clear a mark,
* 'f' - full screen on/off,
* 'q' - quit

The marking labels are saved in the original file along with the synchrotron peak frequency and are available during the following script runnings. One can also change and clean the labels.