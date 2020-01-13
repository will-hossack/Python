import matplotlib.backends
import os.path

be = os.path.dirname(matplotlib.backends.__file__)

print(os.listdir(be))

