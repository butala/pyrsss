import time


class Timer(object):
    def __init__(self):
        self.now = time.time()

    def start(self):
        self.now = time.time()
        return self.now

    def stop(self):
        return time.time() - self.now
