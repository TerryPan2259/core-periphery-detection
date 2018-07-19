from __future__ import absolute_import

import sys
if sys.version_info[:2] < (2, 7):
    m = "Python 3.4 or later is required for NetworkX (%d.%d detected)."
    raise ImportError(m % sys.version_info[:2])
del sys

# Release data
__author__ = '%s <%s>\n%s <%s>\n%s <%s>' % \
    (release.authors['Hagberg'] + release.authors['Schult'] +
        release.authors['Swart'])
__license__ = release.license

__version__ = 0.1 




from .CPAlgorithm import *
from .BE import *
from .MINRES import *
from .KM_config import *
from .KM_ER import *
from .KM_modmat import *
from .qstest import *
from .Cucuringu import * 
from .Rombach import * 
from .Rossa import * 
from .SBM import * 

#from cpalgorithm import *
