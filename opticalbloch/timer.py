import time

class Timer(object):
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print('{:s} done.'.format(self.name))
        print('Time elapsed: {:.3f}s'.format(time.time() - self.tstart))