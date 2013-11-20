
from setuptools import setup, Extension

#############################################################################
### check for Python.h
#############################################################################
from distutils.command.config import config as _config

# a faux config command that checks for the Python.h header file
class config(_config):
  def run(self):
    #self.check_python_dev()
    _config.run(self)

  def check_python_dev(self):
    from distutils import sysconfig
    ok = self.check_header('Python.h',include_dirs=[sysconfig.get_python_inc()])
    if not ok:
      from distutils.errors import DistutilsPlatformError
      errmsg = ("The compiler cannot find the 'Python.h' header file.\n"
                "Please check your configuration to see if you have python-dev installed.")
      raise DistutilsPlatformError(errmsg)

# a faux build_ext command that calls the faux config command
from distutils.command.build_ext import build_ext as _build_ext
class build_ext(_build_ext):
  def run(self):
    _build_ext.run(self)
    print(self.extensions)
    self.build_exec(self.extensions)

  # this is ripped straight from distutils/command/build_clib.py
  def build_exec(self, extensions):
    for ext in extensions:
      exec_name = ext.name
      
      sources = ext.sources #build_info.get('sources')
      if sources is None or not isinstance(sources, (list, tuple)):
        raise DistutilsSetupError(
                "in 'extensions' option (extension '%s'), "
                "'sources' must be present and must be "
                "a list of source filenames" % exec_name)
      sources = list(sources)

      from distutils import log
      log.info("building '%s' objects", exec_name)

      # First, compile the source code to object files in the library
      # directory.  (This should probably change to putting object
      # files in a temporary build directory.)
      extra_compile_args = ext.extra_compile_args or [] #build_info.get('extra_compile_args')

      macros = ext.define_macros[:] # or []build_info.get('macros')
      for undef in ext.undef_macros:
        macros.append((undef,))
      
      include_dirs = ext.include_dirs #build_info.get('include_dirs')
      
      extra_link_args = ext.extra_link_args #build_info.get('extra_link_args')

      objects = self.compiler.compile(sources,
                                      output_dir=self.build_temp,
                                      macros=macros,
                                      include_dirs=include_dirs,
                                      debug=self.debug,
                                      extra_preargs=extra_compile_args,)

      executable = self.compiler.link_executable(objects,
                                                 exec_name,
                                                 output_dir=self.build_temp,
                                                 debug=self.debug,
                                                 extra_postargs=extra_link_args,)

#############################################################################
#### the muritz setup
#############################################################################

muritz = Extension('muritz', 
                   sources = ['/'.join(['src/C',f]) for f in ['alignment.cpp',
                                                              'muritz.cpp',
                                                              'network.cpp',
                                                              'roles.cpp',
                                                              'simulated_annealing.cpp',
                                                             ]
                               ],
                    include_dirs = ['src/C', '/usr/include'],
                    language = 'c++',
                    extra_compile_args = ["-O3",],
                    extra_link_args = ['-lgsl', '-lgslcblas', '-lm', '-lstdc++', ]
                  )

setup(
    name = "muritz",
    version = "0.1",
    description = "muritz is directed graph alignment",
    author = "Daniel B. Stouffer",
    author_email = "daniel.stouffer@canterbury.ac.nz",
    url = 'http://github.com/stoufferlab/muritz',
    packages = ['muritz'],
    package_dir={'muritz':'src/python'},
    ext_modules = [muritz,], # note that this isn't actually a library (as defined by the build_clib hack above)
    cmdclass={'build_ext': build_ext,
              'config': config,
              },
    #scripts=['scripts/muritz'],
    #install_requires=['pymfinder'],
    #dependency_links=["https://github.com/stoufferlab/pymfinder/tarball/master#egg=pymfinder-0.23"]
    #test_suite = 'tests',
)
