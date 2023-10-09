import time

############ To print a counter ############
# Common module used by all files #

# delay in seconds
class Print_erase:
    """
        When the function 'print(s)' is called, it erased the last print and prints 's'
        When the delay between two calls of 'print' is < 'delay', we don't print the second call
        """
    def __init__(self, delay = 1):
        self.last_time = 0
        self.delay = delay

    def __del__(self):
        print()

    def print(self,s):
        if self.delay < (time.time() - self.last_time):
            print("\r", s, end='', flush = True)
            self.last_time = time.time()
