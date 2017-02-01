# Automated documentation using ReadTheDocs and docstrings

The basic documentation for this project is built using the www.readthedocs.org site. It builds documentation from document strings found in the class and function definitions within the LSDPlottingTools package and its submodules.

We use the Google Dosctrings style: http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html

## Some annoying things to watch out for when maintaining the docs

### Module names within LSDPlottingTools
If you change the module names, you must change them in the LSDPlottingTools.rst file, same goes if you add a new module too, it must be either added manually to this file (or you can rebuild with sphinx but I think just adding the module name is easier.)

### External module names
If you have external modules (gdal, cartopy etc) these have to go in the `'MagicMock'` class sys function in `conf.py`. There is a big list of module names that are from external modules that we use in LSDPlottingTools. ReadTheDocs does not have all of these modules installed in the virtualenv it creates so the MagickMock thing tricks it into thinking the modules exist.

The weird thing is that they have to be listed in hierarchical order: e.g. `['gdal', gdal.gdal_array', 'gdal.gdal_array.someothesubmodule']`.

Also weirdly, if you use a 3rd level submodule, but not the second level one e.g. `'cartopy.mpl.patch'` but not `'cartopy.mpl'`, you still have to add `'cartopy.mpl'` to the list and in the correct order (biggest to smallest)

Isn't it so simple and automated...!
