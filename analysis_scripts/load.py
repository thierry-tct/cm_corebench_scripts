
# load data for a particular fault
# Input: fault's res dir tar file
# Output: all tests, tests revealing fault, relevant mutants to tests, mutants to killing test map, tests to killed mutants map

try:
    import muteria
except ImportError:
    sys.path.insert(0, '/home/mfi/mytools/muteria')
    import muteria
    sys.path.pop(0)

def load(tarfile, fault_revealing=True):
    # load fault revealing tests

    # load all tests

    # load execution outputs

    # Compute relevant mutants

    # load matrices and compute mutant killtest mapping

    # remove untared

    # return data
#~ def load()
