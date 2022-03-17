"""
Access stuff from 4FGL
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from pathlib import Path

class PDF(object):
    def __init__(self, pdf_file, size=(500,500)):
        self.pdf = pdf_file
        self.size = size

    def _repr_html_(self):
        return '<iframe src={0} width={1[0]} height={1[1]}></iframe>'.format(self.pdf, self.size)


class SpecPlot(object):
    """Manage access to the 4FGL SED plots
    """
    pattern = '4FGL-DR3_SpecPlots_v*.tgz'
    local_path='sed_plots'

    def __init__(self, cat_path ='/home/burnett/fermi/catalog') -> None:
        import tarfile
        try:
            spec_plots = list(Path(cat_path).glob(self.pattern))[-1] 
        except:
            print(f'Did not find {self.pattern} at {cat_path}')
            return
        self.version = spec_plots.name.split('_')[-1][1:3]
        self.tar = tarfile.open(spec_plots)
        os.makedirs(self.local_path, exist_ok=True)

    def __repr__(self):
        return f'4FGL SED plots at {self.tar.name}'
    
    def __call__(self, source_name, size=(500,500)):
        sn = source_name
        fn = f"SpecPlots_v{self.version}/{sn[6:8]}/"\
             f"{sn.replace(' ','_').replace('.','d').replace('+','p').replace('-','m')}_spec.pdf" 

        try:
            z = self.tar.extractfile(fn).read()
        except KeyError as err:
            print(f'Source {source_name} (file {fn}) not found')
            return
                  
        fnx = Path(self.local_path)/(Path(fn).name) # get name
        with open(fnx, 'wb')as out:
            out.write(z)
        return PDF(fnx,size=size)

class FermiCatalog(object):
    """
    """
    def __init__(self) -> None:
        pass
