---
layout: post
title:  "Build an exe using pyinstaller with GitHub Actions"
date:   2020-10-18 12:06:00
categories: software
tags: [github,python]
thumbnail: https://github.blog/wp-content/uploads/2019/08/DL-V2-LinkedIn_FB.png?fit=1200%2C630
---

## Background

[pyinstaller](https://www.pyinstaller.org/) is used to build standalone executables from Python packages. These are useful particularly if you want to deploy your package as a desktop application for users who can't install Python and all the dependencies on their system. To deploy for Windows users you will be building an exe file. If your package is hosted on github you can use Actions to automate the process of building the executable every time some event happens, like a new release or a commit. All that you need to do is create a `.github/workflow` folder in the root of the repository. Inside this you place one or more .yml files that describe each action.

## Workflow yaml file

This example shows a .yml file describing the pyinstaller build. Some points to note:

* This uses the `workflow_dispatch` event trigger which means it's run manually from the github actions page.
* This assumes there is a `build.spec` file in the root folder that will do the pyinstaller build. This will be specific to your package.
* The pyinstaller action by [Jack McKew](https://github.com/JackMcKew) is used to do the actual build (it uses Wine to emulate windows inside Docker).
* The `with` keyword identifies where the source code is.
* The result is saved as an "artifact" which is placed in the workflow page and can be downloaded as a zip file.

```yml
name: Build GUI exe with Pyinstaller

on:
  workflow_dispatch:
    inputs:
      tags:
        description: 'test build tags'
jobs:
  build:

    runs-on: ubuntu-latest

    steps:
    - uses: actions/checkout@v2

    - name: Package Application
      uses: JackMcKew/pyinstaller-action-windows@main
      with:
        path: .

    - uses: actions/upload-artifact@v2
      with:
        name: my-app
        path: dist/windows
```

## Links

* [pyinstaller](https://www.pyinstaller.org/)
* [Github Actions](https://docs.github.com/en/free-pro-team@latest/actions)
* [PyInstaller-Action-Linux](https://github.com/marketplace/actions/pyinstaller-windows)
