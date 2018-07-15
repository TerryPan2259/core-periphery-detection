from __future__ import absolute_import

import sys
if sys.version_info[:2] < (2, 7):
    m = "Python 2.7 or later is required for NetworkX (%d.%d detected)."
    raise ImportError(m % sys.version_info[:2])
del sys

# Release data
from networkx import release

__author__ = '%s <%s>\n%s <%s>\n%s <%s>' % \
    (release.authors['Hagberg'] + release.authors['Schult'] +
        release.authors['Swart'])
__license__ = release.license

__date__ = release.date
__version__ = release.version




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
