""" Functions that add to the ipython display for navigation within a folder
"""

from . ipynb_docgen import show
from pathlib import Path

__all__ = ['show_folder_index', 'back_to_index', 'show']

here = Path('.')
folder_name = str(here.absolute()).split('/')[-1]

def show_folder_index():
    """ For an "index" notebook
    """
    show(f"""# {folder_name}
    Absolute path: {here.absolute()}

    ## Notebooks:
    """)
    for f in here.glob('*.ipynb'):
        n = f.name.split('.')[0]
        if n=='index': continue
        show(f'* <a href="{n}.ipynb">{n}</a>')
        
def back_to_index():
    """ Insert right-justified link back to index.ipynb
    """
    show(f"""
        <h5 align="right"> <a href="index.ipynb">back to {folder_name} index</a> </h5>
    """)
    
