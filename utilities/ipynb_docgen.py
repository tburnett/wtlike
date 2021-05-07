"""
Extracted from jupydoc to support single-cell output with nbdev
Usage:

from utilities.ipynb.docgen import *

def userdoc():
    '''
    text to show, with {} strings.
    '''
    # (code)
    # 
    return locals()
    
nbdoc(userdoc,)
"""
import sys, os, shutil, string, pprint, datetime
import nbdev # only for a show_doc

__all__ = ['nbdoc', 'image', 'figure', 'monospace', 'capture_print', 'shell', 'create_file'] #,'show_doc']

def doc_formatter(
        text:'text string to process',
        vars:'variable dict'={}, 
        mimetype='text/markdown',
    )->'MimeBundleObject':
    # Returns an object that can be displayed by IPython, interpreted as the mimetype

    # Use a string.Formatter subclass to ignore bracketed names that are not found
    #adapted from  https://stackoverflow.com/questions/3536303/python-string-format-suppress-silent-keyerror-indexerror

    class Formatter(string.Formatter):
        class Unformatted:
            def __init__(self, key):
                self.key = key
            def format(self, format_spec):
                return "{{{}{}}}".format(self.key, ":" + format_spec if format_spec else "")

        def vformat(self, format_string,  kwargs):
            try:
                return super().vformat(format_string, [], kwargs)
            except AttributeError as msg:
                return f'Docstring formatting failed: {msg.args[0]}'
        def get_value(self, key, args, kwargs):
            return kwargs.get(key, Formatter.Unformatted(key))

        def format_field(self, value, format_spec):
            if isinstance(value, Formatter.Unformatted):
                return value.format(format_spec)
            #print(f'\tformatting {value} with spec {format_spec}') #', object of class {eval(value).__class__}')
            return format(value, format_spec)
               
    docx = Formatter().vformat(text+'\n', vars)  if vars else text     
    # enhances this: docx = text.format(**vars)

    class MimeBundleObject(object):
        def _repr_mimebundle_(self, include=None, exclude=None):
            return {mimetype: docx}

    return MimeBundleObject()




def image(filename, 
            caption='', 
            width=None,height=None, 
            image_extensions=['.png', '.jpg', '.gif', '.jpeg'],
            )->'a NBimage object that generates HTML':
    error=''


    #  look in local, 'images' or up one, '../images' 
    for image_path in ['.', 'images', '../images']:
        fullfilename= os.path.join(image_path, filename)
        found = os.path.isfile(fullfilename)
        if found: break 
    if not found:
        error = f'Image file {filename} not found.'
        print(error, file=sys.stderr)

    if not error:
        _, ext = os.path.splitext(filename)
        if not ext in image_extensions:
            error = f'File {filename} not an image? "{ext}" not in {image_extensions}' 
            print(error, file=sys.stderr)

    # class that will be recognized, and FigureWrapper will wrape it.
    class NBimage(object):

        def __init__(self ):
            self.width = width
            self.height = height
            self.caption = caption
            self.fullfilename = fullfilename #????

            if error: 
                self._html = f'<b>{error}</b><br>'
                self.failed = True
     
        def savefig(self, docfilename, **kwargs):
            # call back from formatting -- copy from source 
            shutil.copyfile(self.fullfilename, docfilename )
      
    return NBimage()


try:
    import matplotlib.pyplot as plt
    def figure(fig, caption=None, width=None, height=None):
        """ Set the caption, and possibly width and height attributes of a plt.Figuare object
        """
        assert isinstance(fig, plt.Figure), 'Must be a plt.Figure object'
        fig.caption=caption
        if width: fig.width=width
        if height: fig.height=height
        return fig
except:
    def figure(): pass
    plt=None
try: 
    import pandas as pd
except:
    pd=None


# a dict accumulated here, used to initialze set of wrappers for ObjectReplacer
wrappers = {}
                     
class Wrapper(object):
    """Base class for the replacement classes
    Itself uses the str for the object in question
    """
    def __init__(self, *pars, **kwargs):
        self.obj = pars[0]
        self.vars=pars[1]
        self.indent = kwargs.pop('indent', '5%')
        self.replacer = kwargs.pop('replacer')

    def __repr__(self): return str(self)
    def _repr_html_(self): return str(self)
    def __str__(self):
        text = str(self.obj).replace('\n', '\n<br>')
        return f'<p style="margin-left: {self.indent}"><samp>{text}</samp></p>'



class FigureWrapper(Wrapper): 
    
    def __init__(self, *pars, **kwargs): 
        
        super().__init__(*pars, **kwargs)
        self.indent = kwargs.pop('indent', '5%')

        self.fig = fig = self.obj
        self.__dict__.update(fig.__dict__)
        if getattr(fig, 'failed', False): 
            return

        # from kwargs
        self.folder_name=kwargs.pop('folder_name', 'images')
        self.fig_folders=kwargs.pop('fig_folders', self.replacer.document_folders)

        self.replacer.figure_number += 1
        self.number = self.replacer.figure_number
        self.prefix = self.replacer.figure_prefix
        self.fig_class=kwargs.pop('fig_class', 'nbdoc_image') 

        
        for folder in self.fig_folders:
            t = os.path.join(folder,  self.folder_name)
            os.makedirs(t, exist_ok=True)
            assert os.path.isdir(t), f'{t} not found'
 
    def __str__(self):
        
        if not hasattr(self, '_html') :
        
            # only has to do this once:
            fig=self.fig
            n =self.number
            prefix = self.prefix+'_' if self.prefix else ''

            # the caption, which may be absent.
            caption = getattr(fig,'caption', None)
            if caption is not None:
                caption = f'<b>Figure {n}</b>. ' + getattr(fig,'caption', '').format(**self.vars)
                figcaption = f' <figcaption>{caption}</figcaption>'
            else: figcaption=''

            # assign, or get, a filename
            name =fig.filename if hasattr(fig, 'filename')  else f'{prefix}fig_{n:02d}.png'
            fn = os.path.join(self.folder_name, name )
            browser_fn =fn
            
            # actually save it for the document, perhaps both in the local, and document folders
            # note uses rcParams['savefig.pad_inches']
            for folder in self.fig_folders:
                fig.savefig(os.path.join(folder,fn), bbox_inches='tight') #, pad_inches=0.5)#, **fig_kwargs)
            if plt: plt.close(getattr(fig, 'number', None) )
            img_width = f'width={fig.width}' if hasattr(fig,'width') else ''

            # add the HTML as an attribute, to insert the image, including  caption
            self._html =\
                f'<div class="{self.fig_class}">\n'\
                    f'<figure style="margin-left: {self.indent}" title="Figure {n}">'\
                    f'  <a href="{browser_fn}" title="{browser_fn}">'\
                    f'    <img src="{browser_fn}" alt="Figure {n} at {browser_fn}" {img_width}>'\
                     '   </a>'\
                    f' {figcaption}' \
                     '</figure>\n'\
                '</div>\n'
            #print(f'HTML:\n{self._html}')
        return self._html

    # def __str__(self):
    #     return str(self.img)

if plt: wrappers['Figure']  = (FigureWrapper, {} )
wrappers['NBimage'] = (FigureWrapper, {} )

if pd:
    class DataFrameWrapper(Wrapper): 
        def __init__(self, *pars, **kwargs):

            super().__init__(*pars, **kwargs)
            self._df = self.obj
            kwargs.pop('replacer') # rest should be OK
            self.kw = kwargs

        def __str__(self):
            if not hasattr(self, '_html'):
                self._html = self._df.to_html(**self.kw)                
            return self._html
    df_kwargs= dict( notebook=True, 
                    max_rows=6, 
                    index=False,
                    show_dimensions=False, #True, 
                    justify='right',
                    float_format=lambda x: f'{x:.3f}',
                    )
    wrappers['DataFrame'] = (DataFrameWrapper,  df_kwargs) 

class PPWrapper(Wrapper):
    """Use PrettyPrint
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __str__(self):
        pp = pprint.PrettyPrinter(indent=2)
        text = pp.pformat(self.obj)#.replace('\n', '<br>\n')
        return f'<p style="margin-left: {self.indent}"><samp>{text}</samp></p>'

wrappers['dict'] = (PPWrapper, {} )
wrappers['list'] = (PPWrapper, {} )

class ImageWrapper(Wrapper):
    """ Wrap IPython.display.Image
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._image = self.obj
    def __str__(self):
        return self._image._repr_mimebundle_()
# Placeholder until I figure out how IPython does this in a notebook
# Until then, must create jpeg or png, save with the document, link in in
# wrappers['Image'] = (ImageWrapkper, {})

class ObjectReplacer(dict):
    """
    Functor that will replace objects in a variables dictionary
    It is a dictionary,
        key: a name of a class to have its instances replaced
        value: tuple with two elements: 
            1. the replacement class, which implements a __str__ method
            2. kwargs to apply to new object
    """
    def __init__(self, 
                 folders:'one or more document folders to save images'=['.'], 
                 figure_prefix:'prefix for figure filename'='',
                ):

        self.update(wrappers)
        self.set_folders(folders)
        self.figure_number=0
        self.figure_prefix = figure_prefix
        self.debug=False
        
   
    def set_folders(self, folders):
        # folder management for these guys
        #global document_folders
        self.document_folders = folders
        self.clear()

    def clear(self):
        self.figure_number= 0
 
    @property
    def folders(self):
        return self.document_folders

    def __call__(self, vars):
        """for each value in the vars dict, replace it with a new object that
        implements return of appropriate HTML for the original object
        (Note uses the *class name*, which may not be unique, as a key)
        """
        for key,value in vars.items():
            tkey = value.__class__.__name__
            if self.debug:
                print(f'{key}: {tkey} ')

            new_class, kwargs = self.get(tkey, (None,None))
            
            if new_class:
                newvalue = new_class(value, vars, replacer=self, **kwargs)
                vars[key] = newvalue
                
    def test(self, var:'any object'):
        """Test replacement for a given value. print str(var) before and after replacement
        """
        _class = var.__class__
        _name = _class.__name__
        if _name not in self:
            print(f'no subsitution for class {_name}: "{var}"')
            return

        print(f'replacement: {self[_name]}')
        print(f'{"-"*37}before{"-"*37}\n{var}\n')
        # make a simple vars dict: key is the variable name, value its object
        vars = {'x': var}
        self(vars)
        x = vars['x']
        print(f'{"-"*37}after {"-"*37}\n{x}\n{"-"*80}\n')
        return 


#--------------------External interface------------------------------------------------------------

def monospace(text:'Either a string, or an object',
                summary:'string for <details>'=None,
                open:'initially show details'=False, 
                indent='5%',
                )->str:

    text = str(text).replace('<','&lt;').replace('>','&gt;').replace('\n', '<br>')
    out = f'<div style="margin-left: {indent}"><pre>{text}</pre></div>'
    if not summary:
        return out

    # Set up a "details" HTML tag
    return f'<details {"open" if open else ""} class="nbdoc-description" >'\
           f'  <summary> {summary} </summary>'\
           f'  {out}'\
            ' </details>'
    
def shell(text:'a shell command ', mono=True, **kwargs):
    import subprocess
    try:
        ret = subprocess.check_output([text], shell=True).decode('utf-8')
    except Exception as e:
        ret = f'Command {text} failed : {e}'
    return monospace(ret, **kwargs) if mono else ret

def capture_print(summary=None, **kwargs):
    """
    """

    class Capture_print(object):
        _stream = 'stdout'
        
        def __init__(self):
            import io
            self._new = io.StringIO()
            self._old = getattr(sys, self._stream)

        def __enter__(self):
            setattr(sys, self._stream, self._new)
            return self
        
        def __exit__(self, exctype, excinst, exctb):
            setattr(sys, self._stream, self._old)
            
        def __str__(self):
            return monospace(self._new.getvalue(), summary=summary, **kwargs)

    return Capture_print()

def create_file(func, filename, folder='images'):
    """
    Run func to save filename in images and ../docs/images or docs/images
    return markdown link

    """
    import shutil
    from pathlib import Path
    assert Path('images').is_dir(), 'expect folder {folder} to exist'
    ifilename = Path('images')/filename
    
    try:
        func(ifilename)
    except Exception as msg:
        print(f'Failed to run {func} on {filename}: {msg}')
        raise
            
                
    assert os.path.isfile(ifilename), f'File {filename} was not created in images .'
    ## make sure only goes to folder above this
    #     for p in (Path('../docs'), Path('docs')):
    #         if p.is_dir():
    #             shutil.copyfile(ifilename, p/ifilename)
    shutil.copyfile(ifilename, Path('../docs'/ifilename))
    return f'[{filename}]({ifilename})'    

# convenient interface to show_doc, with disp set to false
def show_doc(elt, **kwargs):
    kwargs.update(disp=False)
    return nbdev.showdoc.show_doc(elt,  **kwargs)

def nbdoc(fun, *pars, name=None, **kwargs):
    """Format the output from an IPython notebook cell using the functon's docstring and computed variables.
     
    If name is specified, use it instead of the function name to distinguish figure file names, say for separate
    executions with differing parameters.

    The required docstring will be interpreted as markdown by IPython.

    args and kwargs will be passed to the user function -- a way to pass information from the notebook environment
    The function must end with "return locals()".
    """
    import inspect
    import IPython.display as display

    # the the docstring and function name
    rawdoc = fun.__doc__
    if rawdoc is None:
        print(f'Function {fun.__name__} must have a docstring', file=sys.stderr)
        return

    doc = inspect.cleandoc(rawdoc)
    name = name or fun.__name__

    # predefine convenient symbols
    vars = dict(date =str(datetime.datetime.now())[:16],
            )

    # run it and collect its local symbol table
    try:
        user_vars = fun(*pars, **kwargs)
    except Exception as e:
        print(f'Function {fun.__name__} failed: {e}', file=sys.stderr)
        raise 
        #return

    if user_vars is None or  type(user_vars)!=dict:
        print( 'The function {fun.__name__} must end with "return locals()"', file=sys.stderr)
        return

    # add convenient symbols
    vars.update(user_vars)

    # check location. Expect the
    if os.path.isdir('docs'):
        # in the root
        folders = ['.', 'docs']
    elif os.path.isdir('../docs'):
        folders = ['.', '../docs']
    else:
        raise Exception('did not find the "docs" folder')
    # initialze the ObjectReplacer
    orep = ObjectReplacer(folders=folders, figure_prefix=name)

    # replace variable objects that are recognized
    orep(vars)

    # format the doc string replacing each recognized "{name}" with a str(obj), where "name" is a key to the object obj
    # in the local symbol table
    md_data = doc_formatter(doc, vars)

    # have IPython display the generated markdown
    display.display( md_data )  
