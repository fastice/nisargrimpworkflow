#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wrap NISAR HDF5 files in a frame directory with VRTs that point directly
into the HDF5 bands, without extracting data to disk.

The VRT files mirror what ``nisarh5toimage --vrtOnly`` produces but are
written via the Python API so no subprocess is spawned.

For RUNW/RIFG the output is a single ``<h5root>.vrt``.
For ROFF each layer gets its own ``<h5root>.layer{N}.vrt`` because
``writeVrtOnly`` in nisarhdf requires a single-layer call and the ROFF
``writeData`` dispatch does not have a ``vrtOnly`` path.
"""

import os
import glob
import nisarhdf
from nisarhdf.writeVrtOnly import writeVrtOnly


# Fields to include in each VRT — mirrors nisarh5toimage 'all fields' lists.
_RUNW_FIELDS = ['unwrappedPhase', 'ionospherePhaseScreen',
                'ionospherePhaseScreenUncertainty', 'coherenceMagnitude',
                'connectedComponents']
_RIFG_FIELDS = ['wrappedInterferogram', 'coherenceMagnitude']
_ROFF_FIELDS = ['slantRangeOffset', 'alongTrackOffset',
                'correlationSurfacePeak', 'snr']
_ROFF_LAYERS = ['layer1', 'layer2', 'layer3']


def wrapH5sInFrameDir(orbit1, frame, verbose=False):
    """
    Wrap all NISAR HDF5 files in ``<orbit1>_<frame>/`` with VRTs that point
    directly into the HDF5 bands.  No data is read into memory.

    Parameters
    ----------
    orbit1 : int or str
        Reference orbit number (used to build the frame directory name).
    frame : int or str
        Frame number.
    verbose : bool, optional
        Print progress messages. The default is False.
    """
    frameDir = f'{orbit1}_{frame}'
    for h5File in sorted(glob.glob(f'{frameDir}/H5/NISAR*.h5')):
        _wrapOneH5(h5File, verbose=verbose)


def _wrapOneH5(h5File, verbose=False):
    """Detect product type and dispatch to the appropriate VRT writer."""
    h5Upper = os.path.basename(h5File).upper()
    for pt in ['RUNW', 'ROFF', 'RIFG']:
        if pt in h5Upper:
            productType = pt
            break
    else:
        if verbose:
            print(f'    wrapH5WithVRT: skipping {h5File} '
                  f'(unrecognised product type)')
        return

    if verbose:
        print(f'    Wrapping {h5File} ({productType}) ...')
    # Use the absolute path so the VRT's HDF5 reference survives chdir
    outputRoot = os.path.splitext(os.path.abspath(h5File))[0]

    if productType == 'RUNW':
        _wrapRUNW(h5File, outputRoot, verbose=verbose)
    elif productType == 'ROFF':
        _wrapROFF(h5File, outputRoot, verbose=verbose)
    elif productType == 'RIFG':
        _wrapRIFG(h5File, outputRoot, verbose=verbose)


def _wrapRUNW(h5File, outputRoot, verbose=False):
    """Open RUNW with noLoadData=True and write a single multi-band VRT."""
    p = nisarhdf.nisarRUNWHDF()
    p.openHDF(h5File, productType='interferogram', noLoadData=True)
    writeVrtOnly(p, outputRoot, _RUNW_FIELDS, verbose=verbose)


def _wrapROFF(h5File, outputRoot, verbose=False):
    """Open ROFF with noLoadData=True and write one VRT per offset layer."""
    p = nisarhdf.nisarROFFHDF()
    p.openHDF(h5File, productType='pixelOffsets', noLoadData=True,
              layers=_ROFF_LAYERS)
    for layer in _ROFF_LAYERS:
        writeVrtOnly(p, f'{outputRoot}.{layer}', _ROFF_FIELDS, layer=layer,
                     verbose=verbose)


def _wrapRIFG(h5File, outputRoot, verbose=False):
    """Open RIFG with noLoadData=True and write a single multi-band VRT."""
    p = nisarhdf.nisarRIFGHDF()
    p.openHDF(h5File, productType='interferogram', noLoadData=True)
    writeVrtOnly(p, outputRoot, _RIFG_FIELDS, verbose=verbose)
