"""
Support single-cell documents in Jupyterlab

Usage: Wtih code like this in a cell

    from utilities.ipynb.docgen import *

    @ipynb_doc
    def userdoc():
        '''
        Markdown text  with {} strings.
        '''
        # (code)
        # 
        return locals()

    userdoc()

one can display the docstring text, with the {} strings replaceed by representations of the variables, like an f-string, but actually equivalent to '...'.format(locals()).
Note the decorator, and the "return locals()" as the last line of the function.


Useful routines:

- image 
- figure
- monospace
- capture_show, capture_hide

- shell
- create_file

Names that are recognized:

- date

About image storage

Images are made from plots, or imported directly, with HTML created to include the image. They are stored in a local folder "images".
The file name by default is the name of the function, or via "userdoc(name='mydoc')". 
"""
import sys, os, shutil, string, pprint, datetime


__all__ = ['nbdoc', 'image', 'figure', 'monospace', 'capture', 'capture_hide', 'capture_show', 'shell', 'create_file', 'ipynb_doc',
        'get_nb_namespace', 'special_prefix', 'figure_number'] #,'show_doc']

special_prefix = ''

# Manage the figure number
class FigureNumber:
    def __init__(self, v=0):
        self.set(v)
    def next(self):
        self.fignum += 1
        return  self.fignum
    def __repr__(self): return f'Current FigureNumber: {self.fignum}'
    def set(self, v=0):
        self.fignum=v
            

    @property
    def value(self):
        return self.fignum

figure_number = FigureNumber(0)

# The decorator to run nbdoc on a function

def ipynb_doc(func):
    def inner(*args, **kwargs):
        nbdoc(func, *args, **kwargs)
    return inner

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
            except Exception as msg:
                return f'Docstring formatting failed: {msg}'
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
    def figure(fig, caption=None, width=None,  base64=True, num=0):
        """ Set the caption, and possibly width of a plt.Figuare object
            width -- in pixels - if set, adjust size for pixels
            num -- set figure number starting value
            base64 -- [True] insert as a "base64" map of pixels. If False, write out as a png file, insert link
        """
        assert isinstance(fig, plt.Figure), 'Must be a plt.Figure object'
        fig.caption=caption
        fig.base64 = base64
        fig.width=width

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
        self.base64 = kwargs.pop('base64', getattr(self.obj, 'base64', True))

        self.fig = fig = self.obj
        self.__dict__.update(fig.__dict__)
        if getattr(fig, 'failed', False): 
            return

        # from kwargs
        self.folder_name=kwargs.pop('folder_name', 'images')
        self.fig_folders=kwargs.pop('fig_folders', self.replacer.document_folders)


        self.prefix = self.replacer.figure_prefix
        self.fig_class=kwargs.pop('fig_class', 'nbdoc_image') 

        
        for folder in self.fig_folders:
            t = os.path.join(folder,  self.folder_name)
            os.makedirs(t, exist_ok=True)
            assert os.path.isdir(t), f'{t} not found'
 
    def __str__(self):
        from IPython.core import pylabtools
        
        if not hasattr(self, '_html') :
        
            # only has to do this once:
            fig=self.fig
            n =  self.number = figure_number.next()
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
            
            if plt: plt.close(getattr(fig, 'number', None) )

            # add the HTML as an attribute, to insert the image, either as base64 or a ref to png file, including  caption

            if self.base64:
                width = getattr(fig,'width', None)
                if  width is not None:
                    # adjust size to match width spec
                    size_inches = fig.get_size_inches()
                    wpix = size_inches[0] * fig.get_dpi(); 
                    fig.set_size_inches(size_inches*fig.width/wpix)
               
                # use IPython tool to create the base64 string for the image
                b64 = pylabtools.print_figure(fig, base64=True, facecolor='white')
                if width is not None:
                    fig.set_size_inches(size_inches)
                self._html =\
                    f'<figure style="margin-left: {self.indent}" title="Figure {n}">'\
                    f'   <img src="data:image/png;base64,{b64}" alt="Figure {n}" '\
                    f' <br> {figcaption}' \
                     '</figure>'
        
            else:
                # actually save it for the document, perhaps both in the local, and document folders
                # note uses rcParams['savefig.pad_inches']
                img_width = f'width={fig.width}' if hasattr(fig,'width') else ''

                for folder in self.fig_folders:
                    fig.savefig(os.path.join(folder,fn), bbox_inches='tight') #, pad_inches=0.5)#, **fig_kwargs)
                self._html =\
                    f'<figure style="margin-left: {self.indent}" title="Figure {n}">'\
                    f'  <a href="{browser_fn}" title="{browser_fn}">'\
                    f'    <img src="{browser_fn}" alt="Figure {n} at {browser_fn}" {img_width}>'\
                     '   </a>'\
                    f' {figcaption}' \
                     '</figure>'
                
        return self._html


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
        self.figure_prefix = figure_prefix
        self.debug=False
        
   
    def set_folders(self, folders):
        # folder management for these guys
        #global document_folders
        self.document_folders = folders
#         self.clear()


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
                show:'initially show details'=False, 
                indent='5%',
                style='',
                )->str:

    text = str(text).replace('<','&lt;').replace('>','&gt;').replace('\n', '<br>')
    out = f'<div style="margin-left: {indent};{style}"><pre>{text}</pre></div>'
    if not summary:
        return out

    # Set up a "details" HTML tag
    return f'<details {"open" if show else ""} class="nbdoc-description" >'\
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

def capture(summary=None, **kwargs):
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
def capture_hide(summary=None, **kwargs):
    return capture(summary, **kwargs)
def capture_show(summary=None, **kwargs):
    return capture(summary, show=True, **kwargs)

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

# # convenient interface to show_doc, with disp set to false
# def show_doc(elt, **kwargs):
#     kwargs.update(disp=False)
#     return nbdev.showdoc.show_doc(elt,  **kwargs)

def get_nb_namespace():
    """
    Return a dictionary, excluding initial underscored variables, of the current IPython notebook namespace
    """
    from IPython.core.interactiveshell import InteractiveShell
    gns = InteractiveShell._instance.get_ipython().user_global_ns
    return {k:v for k,v in filter(lambda t: not t[0].startswith('_'), gns.items())}

def nbdoc(fun, *pars, name=None, fignum=None, **kwargs):
    """Format the output from an IPython notebook cell using the function's docstring and computed variables.
     
    - fun -- User function, which must have a docstring, which will be interpreted as markdown by IPython, and end with `return locals()` 
    - name -- If  specified, use it instead of the function name to distinguish figure file names, say for separate
    executions with differing parameters.
    - fignum -- if specified, apply to the next figure
    
    - *pars, **kwargs -- passed to `fun`
    
    The required docstring 

    """
    import inspect
    import IPython.display as display

    # make sure callable, get name
    if not callable(fun):
        print('nbdoc arg must be callable', file=sys.stderr)
        return

    fun_name =  getattr(fun, '__name__', 'unnamed_function')
    name = name or fun_name
    
    # set the "current" figure number -- will be increment for next figure

    if fignum is not None: 
        figure_number.set( max(fignum-1,0))
    debug = kwargs.pop('debug', False)
    if debug:
        print('Current figure number: ',figure_number)


    # the the docstring and function name
    rawdoc = fun.__doc__
    if rawdoc is None:
        print(f'Function {fun_name} must have a docstring', file=sys.stderr)
        return

    # strip leading spaces from docstring
    doc = inspect.cleandoc(rawdoc)

    # predefine convenient symbols
    vars = dict(date =str(datetime.datetime.now())[:16],
            )

    # run it, adding globals from the IPython session, and collect its local symbol table
    try:
        globals().update(get_nb_namespace())

        user_vars = fun(*pars, **kwargs)
        # user_vars = eval('fun(*pars, **kwargs)', locals().update(get_nb_namespace()) )

    except Exception as e:
        print(f'Function {fun_name} failed: {e}', file=sys.stderr)
        raise 
        #return

    if user_vars is None or type(user_vars)!=dict:
        print( f'The function {fun_name} must end with "return locals()"', file=sys.stderr)
        return

    # add convenient symbols
    vars.update(user_vars)

    # Will put plots images in a local folder, and also docs if found here or in parent
    folders = ['.']
    if os.path.isdir('docs'):      folders.append( 'docs')
    elif os.path.isdir('../docs'): folders.append('../docs')
    else:
        pass
    #    raise Exception('did not find the "docs" folder')
    # initialze the ObjectReplacer
    orep = ObjectReplacer(folders=folders, figure_prefix=special_prefix+name)

    # replace variable objects that are recognized
    orep(vars)

    # format the doc string replacing each recognized "{name}" with a str(obj), where "name" is a key to the object obj
    # in the local symbol table
    md_data = doc_formatter(doc, vars)

    # have IPython display the generated markdown
    display.display( md_data )  
