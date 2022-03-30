---
layout: post
title:  "batchfilerename - A simple utility for batch file renaming"
date:   2022-03-28 15:03:00
categories: python
tags: [python]
thumbnail: /img/batchfilerename_logo.png
---

## Background

<div style="width: 250px; float:right;">
 <img src="/img/batchfilerename_logo.png" width="200px">
</div>

Changing filenames in bulk is time consuming if done manually. If you are changing the filenames according to some pre-determined pattern then it's much easier to automate the process. You can even do this with the `ren` command [in Windows](https://www.makeuseof.com/tag/batch-rename-mass-delete-files-windows/) though it's run from the command line. On linux you can use smart-file-renamer or a host of [command line methods](https://www.makeuseof.com/batch-rename-files-in-linux/). If you want something simple, `batchfilerename` batch file renaming utility written in pure Python. It provides a graphical dialog that let's you find/replace symbols and add text to file names in bulk.

You can install with:

```
pip install batchfilename
```

Install from github:

```
pip install -e git+https://github.com/dmnfarrell/batchfilename.git#egg=batchfilename
```

Install the snap

```
snap install batchfilename
```

## How to use

Run using the command batchfilename. A window with two panes will appear. Select the folder where the files are to be renamed. On the left the files will be listed and on the right a preview of the renamed files is shown (without full path for ease of viewing). You can then select the symbols to find and replace with, which will be applied to all files. Filter the files to be renamed if needed, *.* means all files. **Always use 'preview' first** to check the results before executing as some file name changes might not be reversible. Though you should be able to reverse the last run using the undo button.

Other features:

* The 'occurrences' option allows you to only replace a specific number of instances of a symbol in the name.
* Undo the previous renaming step (assuming you have not quit the program)
* Recursively load a folder

## Links

* [Github page](https://github.com/dmnfarrell/batchfilerename)
* [Snap](https://snapcraft.io/batchfilerename)
* [smart-file-renamer](smart-file-renamer)
