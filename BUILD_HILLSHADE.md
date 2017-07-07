## Fast hillshade module (LSDPlottingTools/fast_hillshade.pyx)

To build the cython hillshade extension, go into the LSDPlottingTools folder and run:

`python setup_cython.py build_ext --inplace`

You will get a fast_hillshade.so file.

Nothing else should be needed to be done and you can import the hillshade module in the usual way, e.g:

`import fast_hillshade as fhill`
