# Mann Kendall Test and Sen's Slope Esitmate on Daily Precipitation Raster Data

The codes are based on the method explained in the following webpages.<br />
* https://www.real-statistics.com/time-series-analysis/time-series-miscellaneous/mann-kendall-test/
* http://www.real-statistics.com/time-series-analysis/time-series-miscellaneous/sens-slope/

This Python code can take daily precipitation <b>raster</b> data as input and return the results of Mann-Kendall Test such as Z-value, P-value and Trend (+1 for positive trend, -1 for negative trend and 0 for no trend) & the results of Sen's slope estimator such as K-value and Slope as Rasters. 

## How to Use
* Make sure that you have [Spyder from Anaconda](https://www.anaconda.com/distribution/) installed.

* Make sure you have these python libraries installed: gdal, numpy, scipy and tqdm. For guide on installing a library using **conda**, visit [this webpage](https://docs.anaconda.com/anaconda/user-guide/tasks/install-packages/).

* In the python code (mann_sen.py), provide daily precipiation data **folder** in variable inpfol, output folder in variable outfol and confidence value in varaible confidence as follows.
```python
inpfol = r"D:\\input"
outfol = r"E:\\output"
confidence = 0.95
```

* Run the code.
