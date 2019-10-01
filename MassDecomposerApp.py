import sys

arg = sys.argv[1]
from os import path
script = path.join(path.dirname(arg),'MassDecomposerScript.py')
if not path.exists(script):
    print('MassDecomposerScript.py not found in data directory!!! Using default.')
    script = path.join(path.dirname(sys.argv[0]),'MassDecomposerScript.py')
else:
    print('using MassDecomposerScript.py in data directory')
import runpy
runpy.run_path(script)