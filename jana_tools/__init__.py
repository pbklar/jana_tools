from .io_utils import RefFileHandler

"""
Tutorial reference:

package/
  models/
    api.py
    db.py
    user.py
    
from package.models.api import MyAPI
from package.models.db import UserTable    

__init__.py
from package.models.api import MyAPI
from package.models.db import UserTable

Result:
from package import MyAPI, UserTable, ...
"""

