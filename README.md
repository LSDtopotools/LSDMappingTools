# LSDMapArtist
An experimental python package for plotting LSDTopoTools output and world domination. (But mainly plotting)


An object-oriented plotting module for constructing
drape maps in a reusable, generic way.

Experimental. Use at your own risk.

If it ever works I'll transfer it to @LSDTopoTools

This software is realsed under the Artistic Licence v2.0

## Example use:
(subject to change and open to suggestions)

```
import lsdmapartist as lsdmap

myDrapePlot = lsdmap.DrapePlot(RasterName, BackgroundRasterName,
                                Directory, "Blues", drape_min_threshold=0.05)
                                
myDrapePlot.set_drape_min_threshold(0.1)
myDrapePlot.set_colourmap("jet")
myDrapePlot.set_backgroundtype("Hillsahde")
myDrapePlot.set_backgroundtype("TerrainColours")

myDrapePlot.save_fig("png", "ESurfDyn")

```
And so on...
