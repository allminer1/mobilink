"""Project-local alias, not part of the released package.

The package was called linkstat while the analyses of this project were run (frozen commit 9819182 and earlier) and is
released as mobilink. The analysis scripts in scripts/ and frozen/ import `linkstat`; this alias lets them run
unchanged on the renamed package. It is on the path of envs/linkstat_dev through site-packages/linkstat_compat.pth.
"""
import importlib
import pkgutil
import sys

import mobilink as _m

sys.modules[__name__] = _m
for _mod in pkgutil.walk_packages(_m.__path__, "mobilink."):
    sys.modules["linkstat." + _mod.name[len("mobilink."):]] = importlib.import_module(_mod.name)
