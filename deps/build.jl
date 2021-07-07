import PyCall: pyimport

# See https://stackoverflow.com/questions/12332975/installing-python-module-within-code.
const PIP_PACKAGES = ["f90nml", "scipy", "matplotlib"]

sys = pyimport("sys")
subprocess = pyimport("subprocess")
subprocess.check_call([sys.executable, "-m", "pip", "install", "--user", "--upgrade", "--force-reinstall", PIP_PACKAGES...])