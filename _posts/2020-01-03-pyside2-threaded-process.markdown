---
layout: post
title:  "Concurrent processes in PySide2/PyQt5 applications"
date:   2020-01-03 10:04:00
categories: python
tags: [pyside2,python]
thumbnail: /img/pyside2-threading.png
---

## Background

<div style="width: 350px; float:right;">
 <a href="/img/pyside2-threading-example.gif"> <img src="/img/pyside2-threading-example.gif" width="300px"></a>
</div>

Qt is a popular GUI toolkit for writing desktop applications. It has bindings for Python using either PySide2 or PyQt5 (which use essentially identical syntax). This means you can code Qt apps without needing to know C++. One of the challenging aspects in all GUI tookits is running processes concurrently while still using your application. This is because the main GUI runs in a thread and if you launch your process in the same thread it will block all user interaction until it's finished. That's a problem if your task runs more than a second or two. The solution is to run your jobs in other threads. The recommended solution in Qt is to use classes called `QRunnable` and `QThreadPool`. `QThreadPool` handles queuing and execution of workers. Workers are `QRunnable` objects containing the process you want to run. You place the code to run in the `run()` method of the `QRunnable` class.

This method is mostly taken from the [example](https://www.learnpyqt.com/courses/concurrent-execution/multithreading-pyqt-applications-qthreadpool/) by Martin Fitzpatrick at learnpyqt.com. It's also explained better and in more detail there.

## Imports

Imports for PySide2 are somewhat different than the PyQt5 version.

```python
import sys,os,subprocess,time,traceback
import random
from PySide2 import QtCore
from PySide2.QtWidgets import *
from PySide2.QtGui import *
```

## Worker and ThreadPool

 We make a Worker by sub-classing `QRunnable`, then placing the code wewish you execute within the `run()` method. It also defines `Signals` that are used to pass data out of the running Worker. So you can emit progress changes to the GUI. You probably shouldn't use this to pass very large amounts of data back to the GUI though.

```python
class Worker(QtCore.QRunnable):
    """Worker thread for running background tasks."""

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()
        self.kwargs['progress_callback'] = self.signals.progress

    @QtCore.Slot()
    def run(self):
        try:
            result = self.fn(
                *self.args, **self.kwargs,
            )
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()

class WorkerSignals(QtCore.QObject):
    """
    Defines the signals available from a running worker thread.
    Supported signals are:
    finished
        No data
    error
        `tuple` (exctype, value, traceback.format_exc() )
    result
        `object` data returned from processing, anything
    """
    finished = QtCore.Signal()
    error = QtCore.Signal(tuple)
    result = QtCore.Signal(object)
    progress = QtCore.Signal(int)
```

## The application

This is a small dialog to illustrate the principle. We make a QDialog window and add start and stop buttons and a progress bar. The `run_threaded_process()` method defined here can run any process in a thread. It instantiates a Worker and then adds it to the threadpool. `progress_fn()` and `on_complete()` are provided as arguments and connected to the signals of the worker. Note that to stop the process I just added a flag that is set when the stop button is pressed. Finally, `test()` is the actual function we run. It could be anything but in this case it just runs a loop calculating randon numbers. The loop inside `test()` is then interrupted. This was the only way I knew how to do it but it may not be the best method. Also notice that if you keep pressing the start button new processes will run and more numbers appear. The stop button stops all of them at once.

```python
class App(QDialog):
    """GUI Application using PySide2 widgets"""
    def __init__(self):
        QDialog.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setGeometry(QtCore.QRect(200, 200, 500, 500))
        self.threadpool = QtCore.QThreadPool()

        layout = QVBoxLayout(self)
        self.setLayout(layout)
        self.startbutton = QPushButton('START')
        self.startbutton.clicked.connect(self.run)
        layout.addWidget(self.startbutton)
        self.stopbutton = QPushButton('STOP')
        self.stopbutton.clicked.connect(self.stop)
        layout.addWidget(self.stopbutton)
        self.progressbar = QProgressBar(self)
        self.progressbar.setRange(0,1)
        layout.addWidget(self.progressbar)
        self.info = QTextEdit(self)
        self.info.append('Hello')
        layout.addWidget(self.info)
        return

    def progress_fn(self, msg):
        """Update progress"""

        self.info.append(str(msg))        
        return

    def run_threaded_process(self, process, progress_fn, on_complete):
        """Execute a function in the background with a worker"""

        worker = Worker(fn=process)
        self.threadpool.start(worker)
        worker.signals.finished.connect(on_complete)
        worker.signals.progress.connect(progress_fn)
        self.progressbar.setRange(0,0)
        return

    def run(self):
        """call process"""

        self.stopped = False
        self.run_threaded_process(self.test, self.progress_fn, self.completed)

    def stop(self):
        self.stopped=True
        return

    def completed(self):
        self.progressbar.setRange(0,1)
        return

    def test(self, progress_callback):
        """Do some process here"""

        total = 500
        for i in range(0,total):
            time.sleep(.2)
            x = random.randint(1,1e4)
            progress_callback.emit(x)
            if self.stopped == True:
                return
```

## Links

* [Script with this code](https://github.com/dmnfarrell/teaching/blob/master/gui/pyside_threading.py)
* [Multithreading pyqt applications](https://www.learnpyqt.com/courses/concurrent-execution/multithreading-pyqt-applications-qthreadpool/)
