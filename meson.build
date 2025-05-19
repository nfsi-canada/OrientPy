project('orientpy', 'c',
  version : '0.2.1',
  license: 'MIT',
  meson_version: '>=0.64.0',
)

py_mod = import('python')
py = py_mod.find_installation()
py_dep = py.dependency()

py.install_sources(
  'orientpy/__init__.py',
  'orientpy/classes.py',
  'orientpy/plotting.py',
  'orientpy/utils.py',
  'orientpy/io.py',
  'orientpy/scripts/bng_calc.py',
  'orientpy/scripts/bng_average.py',
  'orientpy/scripts/dl_calc.py',
  'orientpy/scripts/dl_average.py',
  'orientpy/scripts/__init__.py',
  pure: false,
  preserve_path: true,
  )
  
  install_data(
    'orientpy/dispmaps/R.gv.10.txt',
    'orientpy/dispmaps/R.gv.15.txt',
    'orientpy/dispmaps/R.gv.20.txt',
    'orientpy/dispmaps/R.gv.25.txt',
    'orientpy/dispmaps/R.gv.30.txt',
    'orientpy/dispmaps/R.gv.35.txt',
    'orientpy/dispmaps/R.gv.40.txt',
    install_dir: py.get_install_dir(pure: false) / 'orientpy/dispmaps',
)
