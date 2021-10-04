"""
Copied to wtlike/config
"""
import time

class Timer():
    """Usage:
    ```
    with Timer() as t:
        time.sleep(5)
    print(t, t.elapsed)
    ```
    """
    def __init__(self):
        self.t=time.time()
        self.exit_time=0
        
    def __enter__(self):
        return self
    def __exit__(self, *pars):
        #print(f'Elapsed time: {self():.0f} sec')
        self.exit_time = time.time()-self.t
    def __repr__(self):
        return f'Elapsed time: {self.elapsed:.1f} sec'
    @property
    def elapsed(self):
        return min(time.time()-self.t, self.exit_time) 
