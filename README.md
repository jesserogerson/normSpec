#normalizeSpectra.py
Author: Jesse A. Rogerson, email: rogerson@yorku.ca

The code implementation was written by the above, but with helpful (and crucial) contributions from: **Patrick B. Hall, Paola Rodriguez Hidalgo**

-----
### Synopsis

These code takes a list of raw spectra from the command line and allows the user
to interactively normalize them. The program was specifically desinged to be
used for the main author's PhD thesis, and thus isn't necessarily general to
all raw spectra.

The program can broken down into roughly two separate functions:

1. normalize()
This function allows the user to interactively set the parameters of
normalization and then execute the normalization.

2. plotNorm()
This method is called automatically after the users has asked these code to
calculate the normalized spectra. plotNorm() allows the user to interactively
augment a plot of the normalized spectra in order to both check if the
normalization parameters were good/not, as well as have a publication-worthy
figure if desired.

After the users elects to quit plotNorm(), they are sent back to the original
normalize() function where they may augment the parameters and redo the
normalization, or quit normalizeSpectra.py altogether.

###To Run:
----------

$> ./normalizeSpectra.py JHHMMSS.card

The *.card file is both the list of raw spectra and where the major Information
of the object is held. In order for these code to run the *.card file must
be structured in the following way:

###The *.card file:
----------
JHHMMSS.card <--- MUST be called this
contents:
line0: SDSS Jhhmmss.ss+/-ddmmss.s
line1: RA Dec
line2: redshift'
line3: label1 MJD1 /path/to/raw/spec1/
line4: label1 MJD2 /path/to/raw/spec2/
line5: label1 MJD3 /path/to/raw/spec3/
line6: label1 MJD4 /path/to/raw/spec4/
                .
lineF: labelN MJDN /path/to/raw/specN/

The *.card file may have an arbitrary number of raw spectra for normalization
but the first three lines MUST be name/RA DEC/redshift. The raw spectra lines
MUST be label MHD path (space separated).


### Installation

No Installation needed, just download and execute script.

### License

see LICENSE.txt
