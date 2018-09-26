---
layout: post
title:  "Deploy a Python application with snapcraft"
date:   2018-09-26 10:12:00
categories: software
tags: [snapcraft]
---

## Snaps

Snaps are a convenient way to build binaries for many Linux variants at once. These were introduced by Ubuntu as a way of partly solving the application distribution problem on Linux. Linux has multiple packaging systems used in different versions. This makes it hard for a developer to provide a single simple installer for their program that is platform independent (one of the few advantages of using Windows). Snaps (and also flatpak) work by packaging all the dependencies needed together and confining the app to it's own environment, thus also improving security. Once you build your application as a snap you can roll out updates when needed, independent of any Linux distribution. Snaps are supported on most major versions of Linux. For more information about snaps see the links at the end.

## Building a Python snap

This article is about turning your Python application (that has a command line interface or desktop application) into a snap.
Python has it's own tool for packaging called pip which works quite well. But it can still be system dependent. Snapcraft builds on top of this to create a self contained package that you can deploy. This should work consistently as it uses it's own version of python and dependent libraries. You can follow the simplified example below or go to the more detailed [official documentation](https://docs.snapcraft.io/build-snaps/python).

## Simple example

The snapcraft program is used to build snaps. This should be installed on your system. It can be added in Ubuntu using apt:

```sudo apt install snapcraft```

Say you have package called testapp that you want to snap. It would typically have the file structure below.  Like most python packages it has a setup.py file. This is needed to make the snap. You can download these sample files <a href="/other/python-snap-template.zip">here</a> and try it for yourself.

```
├── setup.py
├── snap
│   ├── gui
│   │   └── icon.png
│   ├── snapcraft.yaml
└── src
    ├── app.py
    ├── description.txt
    └── __init__.py
```

Notice that we have created a folder called snap which contains a file *snapcraft.yaml*. This is where you configure your snap build settings.

```
name: testapp
version: '0.1'
#if using a git repository you can use the current revision number
version-script: git describe --abbrev=1 --tags
summary: python test package
description: |
 this is a test snap for python

grade: stable
confinement: devmode

apps:
  testapp:
    command: bin/hello
    plugs: [home,network-bind]   
parts:
  testapp:
    plugin: python
    python-version: python3
    source: ../
    stage-packages: [tcsh,gawk]
    python-packages: [numpy]
```

Briefly, each .yaml file has a section with the name of the snap, a summary and description. Then there are *parts* which can be anything you want added such as programs, libraries, or other assets needed. The above example uses one Python part to fetch the source from the parent directory. This could also be a github repository or URL. This uses the setup.py file to build the required python packages. You can optionally add a stage-packages line that lists other programs to be added from the apt repository. The python-packages line is only needed if you want extra python packages not in the setup.py file.
The *apps* section allows you to expose commands that will be available when the snap is installed. In this case testapp is a command (does not have to be the same as the snap name) is a link to bin/hello which is exposed in the setup.py file as shown below:

```
entry_points = {
    'console_scripts': [
        'hello=src.app:main']
        },
```

If your command name matches the snap name, users will be able run the command directly.

The [confinement](https://docs.snapcraft.io/reference/confinement) mode determines how self contained the snap is. Generally strict should be used. You can used devmode as in this case to test things until you have it working.

To build the snap, change to the snap directory and type:

```snapcraft snap```

## Deploying

When you build a file like testapp_latest_amd64.snap will be created. This is the final snap. It can in principle be given to anyone else to install locally. If you distribute that way the user will have to add the --dangerous option when installing which simply means the file has not been verfied via the store. If you want many people to see and use the snap you really need to use the store.

Ubuntu provides the [snap store](https://snapcraft.io/store) infrastructure making it very easy to deploy and host your package. Though you do have to have an account on ubuntu one to use it.

When the snap is ready for putting on the store you can follow these steps:

1. change the confinement to strict and rebuild the snap
2. publish using these steps:

```
snapcraft login #requires you are registered on ubuntu one
snapcraft register testapp # done once only
snapcraft push testapp_0.1_amd64.snap
snapcraft release testapp 1 stable #where 1 is the revision number
```

Once this is done you can check the snap store for the application. It will provide instructions for users to install it themselves. This is a brief overview. More detail on each step is available in the documentation.

## Links

* [Python snaps](https://docs.snapcraft.io/build-snaps/python)
* [Canonical’s Snap: The Good, the Bad and the Ugly](https://thenewstack.io/canonicals-snap-great-good-bad-ugly/)
* [snapcraft home](https://snapcraft.io/build)
* [snapcraft forum](https://forum.snapcraft.io/)
* [flatpak](https://flatpak.org/)
