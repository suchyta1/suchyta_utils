"""
The :mod:`system` submodule has functions for working with system resources.
So far it's only memory usage, but I might add more.

"""

import numpy as _np
import subprocess as _subprocess
import os as _os
import resource as _resource


def _GetDen(unit):
    power = 2
    if unit.upper() in ['M','MB']:
        power = 1
    elif unit.upper() in ['K','KB']:
        power = 0
    den = _np.power(1024.0, power)
    return den


def GetMaxMemoryUsage(unit='GB'):
    """
    Check the maximum memory usage by the current python process.
    This function uses the python :mod:`resource` module.

    Parameters
    ----------
    unit (str)
        The unit to return the memory usage: 
            | 'G'/'GB' for gigabytes,
            | 'M'/'MB' for megabytes, 
            | 'K'/'KB' for kilobytes

    Returns
    -------
    memory (float)
        Max memory usage in the chosen units

    """

    den = _GetDen(unit)
    mem = _resource.getrusage(_resource.RUSAGE_SELF).ru_maxrss / den
    return mem


def GetCurrentMemoryUsage(unit='GB'):
    """
    Check how much memory the current python process is using.
    This function uses the unix :mod:`ps` command.

    Parameters
    ----------
    unit (str)
        The unit to return the memory usage: 
            | 'G'/'GB' for gigabytes,
            | 'M'/'MB' for megabytes, 
            | 'K'/'KB' for kilobytes

    Returns
    -------
    memory (float)
        Memory usage in the chosen units

    """

    den = _GetDen(unit)
    p = _subprocess.Popen(['ps', 'v', '-p', str(_os.getpid())],stdout=_subprocess.PIPE)
    out = p.communicate()[0].split(b'\n')
    vsz_index = out[0].split().index(b'RSS')
    mem = float(out[1].split()[vsz_index]) / den
    return mem
