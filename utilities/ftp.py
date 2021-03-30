"""
"""
import os, sys
import shutil
import numpy as np



import pysftp as sftp

class SLAC(object):
       
    server =     'rhel6-64.slac.stanford.edu'
    username =   'burnett'

    def __init__(self,  
            remote_path='', 
            local_path='/tmp',):
        self.local_path = local_path
        self.remote_path = remote_path
        os.makedirs(self.local_path, exist_ok=True)
        self.svr=None

    def get(self, files:'files to download', 
                reload=False, quiet=True):
        """
        supports final wild card: files can end with *
        will run makedirs on files
        """

        not_connected =  self.svr is None 
        if not_connected:    self.connect()

        if type(files) == str and files[-1]=='*':
            # process wildcard
            found = self.svr.listdir(files[:-1])
            # add path
            p,_ = os.path.split(files)
            files = [os.path.join(p,f) for f in found]

        for f in np.atleast_1d(files): 
            _,fext = os.path.splitext(f)
            if not fext: # should check directly for a folder
                continue
            p,q = os.path.split(f)
            if p:
                lp =  os.path.join(self.local_path,p)
                if not os.path.exists(lp):
                    os.makedirs(lp, exist_ok=True)
            locf = os.path.join(self.local_path, f)
            if reload or not os.path.isfile(locf):  
                try:            
                    self.svr.get(f, locf)
                    if not quiet: print(f'ftp: {f} -> {locf}')
                except Exception as msg:
                    print(f'Failed to get {f} to {locf}: {msg}', file=sys.stderr)
                    raise
        
        if not_connected: self.close()

            
    def put(self, files, replace=False, quiet=True):
        not_connected =  self.svr is None 
        if not_connected:    self.connect()

        for f in np.atleast_1d(files):
            locf = os.path.join(self.local_path, f)
            assert os.path.isfile(locf), f' {locf} does not exist'
            try:
                self.svr.put(locf, f) 
                if not quiet: print(f'ftp: {locf} -> {f}')
            except Exception as msg:
                print(f'ftp put, {locf}->{f} failed: {msg}', file=sys.stderr)
                raise

        if not_connected: self.close()

    def chdir(self, path):
        if self.svr is None:
            print(f'Not connected')
            return
        self.svr.chdir(path)

    def connect(self):
        try:
            self.svr =  sftp.Connection(self.server, self.username)
         
        except:
            print(f'Failed to connect to {self.username}@{self.server}',
                file=sys.stderr)
            raise
        try:
            if self.remote_path: 
                self.svr.chdir(self.remote_path)
        except Exception as msg:
            print(f'chdir failed to {self.remote_path}: {msg}', file=sys.stderr)
            raise

    def listdir(self, remote_path=''):
        if self.svr is None:
            print(f'Not connected')
            return
        try:
            return self.svr.listdir(remote_path)
        except Exception as msg:
            print(f'Failded listdir({remote_path}) in {self.remote_path}: {msg}', file= sys.stderr)

    def close(self): 
        if self.svr: self.svr.close()
        self.svr=None

    def __enter__(self):
        self.connect() 
        return self
  
    def __exit__(self,  type, value, tb): 
        self.close() 

    def __str__(self):
        return f'FTP to SLAC: {self.remote_path},\n local: {self.local_path}'\
                f' \n {"NOT" if self.svr is None else ""} connected'
    def __repr__(self): return str(self)

    