# MIET-GUI

MIET_GUI - A graphical user interface to explore MIET using MatLab 

## Description

The discovery of Förster resonance energy transfer (FRET) has revolutionized our ability to measure inter- and intramolecular distances on the nanometre scale using fluorescence imaging. The phenomenon is based on electromagnetic-field-mediated energy transfer from an optically excited donor to an acceptor. We replace the acceptor molecule with a metallic film and use the measured energy transfer efficiency from donor molecules to metal surface plasmons to accurately deduce the distance between the molecules and metal. Like FRET, this makes it possible to localize emitters with nanometre accuracy, but the distance range over which efficient energy transfer takes place is an order of magnitude larger than for conventional FRET. This creates a new way to localize fluorescent entities on a molecular scale, over a distance range of more than 100 nm.

The MIET-GUI (metal-induced energy transfer graphical user interface) is a software packet written in the MATLAB programming language. It is intended as a user-friendly way to convert lifetime-images of arbitrary dipole emitters into height-images. The theory behind this software packet is described in: N. Karedla, D. Ruhlandt, A. M. Chizhik, J. Enderlein, A. I. Chizhik, "Metal-Induced Energy Transfer", Advanced Photon Counting: Applications, Methods, Instrumentation; Springer Series on Fluorescence, edited by P. Kapusta et al.; pp. 265-281 (Springer, 2014).

Briefly, a dipole emitter can transfer energy to nearby metal surfaces, thus changing its excited state lifetime. This near-field effect is distance dependent with an almost linear regime in the first tens of nanometers, depending on the wavelength of the emitter and the properties of the sample. Our software uses the relationship between lifetime and height of the emitter above the surface to transform lifetime information into height information.

## Authors and acknowledgment
This software package includes valuable contributions from 
Daja Ruhland
Alexej Chizhik
Narain Karedla
Sebasitan Isbaner
Jörg Enderlein

## License
For open source projects, say how it is licensed.

