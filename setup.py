
import os
from setuptools import setup, Extension, find_packages
from distutils import log
from distutils.errors import DistutilsPlatformError, DistutilsSetupError
from stat import ST_MODE

#############################################################################
### check for Python.h
#############################################################################
from distutils.command.config import config as _config

# a faux config command that checks for the Python.h header file
class config(_config):
  def run(self):
    self.check_python_dev()
    _config.run(self)

  def check_python_dev(self):
    from distutils.sysconfig import get_python_inc
    log.info("cehcking for 'Python.h' header file")
    ok = self.check_header('Python.h',include_dirs=[get_python_inc()])
    if not ok:
      errmsg = ("The compiler cannot find the 'Python.h' header file.\n"
                "Please check your configuration to see if you have python-dev installed.")
      raise DistutilsPlatformError(errmsg)

# a faux install command since the standard distutils version is missing build_clib
from distutils.command.install_lib import install_lib as _install_lib
class install_lib(_install_lib):
  def run(self):
    self.build()
    outfiles = self.copy_tree(self.build_dir, self.install_dir)
    if outfiles is not None and self.distribution.has_pure_modules():
      self.byte_compile(outfiles)

  def build(self):
    if not self.skip_build:
      if self.distribution.has_pure_modules():
        self.run_command('build_py')
      if self.distribution.has_ext_modules():
        self.run_command('build_py')
      if self.distribution.has_c_libraries():
        self.run_command('build_clib')

# a faux build_ext command that calls the faux config command
from distutils.command.build_clib import build_clib as _build_clib
class build_clib(_build_clib):
  def run(self):
    if not self.libraries:
      return

    from distutils.ccompiler import new_compiler
    self.compiler = new_compiler(compiler=self.compiler,
                                 dry_run=self.dry_run,
                                 force=self.force)
    from distutils.sysconfig import customize_compiler
    customize_compiler(self.compiler)

    if self.include_dirs is not None:
      self.compiler.set_include_dirs(self.include_dirs)
    if self.define is not None:
      # 'define' option is a list of (name,value) tuples
      for (name,value) in self.define:
        self.compiler.define_macro(name, value)
    if self.undef is not None:
      for macro in self.undef:
        self.compiler.undefine_macro(macro)

    self.build_exec(self.libraries)

  # this is ripped straight from distutils/command/build_clib.py
  def build_exec(self, libraries):
    for (exec_name, build_info) in libraries:

      sources = build_info.get('sources')
      if sources is None or not isinstance(sources, (list, tuple)):
        raise DistutilsSetupError(
                "in 'extensions' option (extension '%s'), "
                "'sources' must be present and must be "
                "a list of source filenames" % exec_name)
      sources = list(sources)

      log.info("building '%s' objects", exec_name)

      macros = build_info.get('macros')
      extra_compile_args = build_info.get('extra_compile_args')
      include_dirs = build_info.get('include_dirs')

      objects = self.compiler.compile(sources,
                                      output_dir=self.build_temp,
                                      macros=macros,
                                      include_dirs=include_dirs,
                                      debug=self.debug,
                                      extra_preargs=extra_compile_args,)


      # DBS: added this to mimic creating a library (but instead create and executable with the same name)
      log.info("building '%s' executable", exec_name)

      extra_link_args = build_info.get('extra_link_args')

      executable = self.compiler.link_executable(objects,
                                                 exec_name + '.x',
                                                 output_dir='./scripts/',
                                                 debug=self.debug,
                                                 extra_postargs=extra_link_args,)

      #from distutils.file_util import copy_file
      #ext_path = self.get_ext_fullpath(ext.name)
      #copy_file('/'.join([self.build_temp,exec_name]), ext_path)

#############################################################################
#### the muritz setup
#############################################################################
# the c++ extension module
extension_mod = Extension("muritzex", 
    sources = ['/'.join(['src/C',f]) for f in ['alignment.cpp',
    'muritz.cpp',
    'network.cpp',
    'roles.cpp',
    'simulated_annealing.cpp',
    'muritz_module.cpp',
    'anneal.cpp']
    ],
    include_dirs = ['src/C', '/usr/include', '/usr/include/python2.7'],
    library_dirs = ['/usr/lib'],
    libraries = ['gsl', 'gslcblas', 'm'], 
    language = 'c++',
    extra_compile_args = ["-O3", "-std=gnu++11"],
    extra_link_args = ['-lgsl', '-lgslcblas', '-lm', '-lstdc++'],
)

setup(
    name = "muritzex",
    ext_modules=[extension_mod]
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
    #libraries = [muritz,], # note that this isn't actually a library (as defined by the build_clib hack above)
    cmdclass={'build_clib': build_clib,
              'config': config,
              'install_lib': install_lib,
              },
    scripts=['scripts/muritz'],
    #install_requires=['pymfinder'],
    #dependency_links=["https://github.com/stoufferlab/pymfinder/tarball/master#egg=pymfinder-0.23"]
    #test_suite = 'tests',
)


