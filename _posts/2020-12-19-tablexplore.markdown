---
layout: post
title:  "Tablexplore - a desktop tool for table analysis"
date:   2020-12-19 12:28:00
categories: software
tags: [PySide2,tablexplore]
thumbnail: /img/tablexplore_logo.png
---

## Background

<div style="width: 240px; float:right;">
<img src="/img/tablexplore_logo.png" width="200px">
</div>

[Tablexplore](https://github.com/dmnfarrell/tablexplore) is an application for data analysis and plotting built in Python using the PySide2/Qt toolkit. The interface allows quick visualization of data with convenient plotting. The primary goal is to let users explore their tables interactively without any prior programming knowledge and make interesting plots as they do this. The program is intended mainly for educational/scientific use but should be useful for a variety of general applications.

This isn't a replacement for a spreadsheet as such. Instead it's meant for quick and easy exploration of tabular data. It can do some things much more conveniently and others not as well. One advantage is the ability to load and work with relatively large tables as compared to spreadsheets. The focus is on data manipulation rather than data entry, though basic cell editing and row/column changes are supported. The program is free and open source. This is an updated version of an older package called [dataexplore](https://github.com/dmnfarrell/pandastable).

### Current features

* spreadsheet-like interface for table manipulation - add/remove rows and columns
* copy and paste between tables
* table analysis tools such as groupby-combine, pivot, merge, join and concatenate
* basic table formatting such as font, cell color, text size and column width
* import/export of supported text files
* rendering of large tables is possible
* simple undo mechanism
* moderately advanced plotting with a simple interface
* data cleaning
* python interpreter to manipulate table if needed or for learning Python/Pandas

## Interface

<div style="width: auto;">
 <a href="/img/tablexplore_scr1.png"> <img class="scaled" src="/img/tablexplore_scr1.png"></a>
  <p class="caption">Main window showing table, plot view and python interpreter. </p>
</div>

<div style="width: auto;">
 <a href="/img/tablexplore_gui.gif"> <img class="small-scaled" src="/img/tablexplore_gui.gif"></a>
  <p class="caption">Brief example of interface usage.</p>
</div>

## Installing

This tool will run on Windows and Linux and probably OS X (not tested). On Linux you can install with pip as follows:

```bash
pip install -e git+https://github.com/dmnfarrell/tablexplore.git#egg=tablexplore
```

There is a standalone binary file for Windows users that can be downloaded [here](https://dmnfarrell.github.io/tablexplore/).

## Using in Python

You can also use the `DataFrameWidget` widget in your own python Qt application. This is shown in the code below.

```python
from PySide2 import QtCore
from PySide2.QtWidgets import *
from PySide2.QtGui import *
import pandas as pd
from tablexplore import data, core, plotting, interpreter

class TestApp(QMainWindow):
    def __init__(self, project_file=None, csv_file=None):

        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("Example")
        self.setGeometry(QtCore.QRect(200, 200, 800, 600))
        self.main = QWidget()
        self.setCentralWidget(self.main)
        layout = QVBoxLayout(self.main)
        df = data.getSampleData()
        t = core.DataFrameWidget(self.main,dataframe=df)
        layout.addWidget(t)
        #show a Python interpreter
        t.showInterpreter()
        return

if __name__ == '__main__':
    import sys
    app = QApplication(sys.argv)
    aw = TestApp()
    aw.show()
    app.exec_()
```

## Technical

This tool uses the [pandas](https://pandas.pydata.org/) DataFrame class to store the table data. Pandas is an open source Python library providing high-performance data structures and data analysis tools.

## Links

* [github page](https://github.com/dmnfarrell/tablexplore)
* [homepage](https://dmnfarrell.github.io/tablexplore/)
* [Thoughts on the R Tidyverse](https://towardsdatascience.com/a-thousand-gadgets-my-thoughts-on-the-r-tidyverse-2441d8504433)
