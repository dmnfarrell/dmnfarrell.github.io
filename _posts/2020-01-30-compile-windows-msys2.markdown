---
layout: post
title:  "Compile windows exe files with MSYS2"
date:   2020-01-30 12:15:00
categories: software
tags: [windows]
thumbnail: /img/msys2.png
---

## Background

Many programs for scientific use are written in the C/C++ languages. They're also developed and run on POSIX based systems such as Linux. They can be cross-compiled for windows systems but those that require Unix/POSIX compliance won't work unless the author has also tailored them for windows use. So what if there is a unix-based program you want to use in windows but it's not available? MSYS2 is a system that can help solve this problem by providing a unix like environment inside windows on which you can compile a program. It provides lots of programs ported from unix like a set of GCC compilers for C and C++. It features a package management system called Pacman to provide easy installation of packages.

There is only one catch to making the .exe this way: it needs a file called **msys-2.0.dll** to be present for it to run. This library provides the features of POSIX that basically make the program think it's not running on windows. If your program doesn't need that kind of support you can use the MINGW-64 shell to compile instead. Then your executable won't need the dll file. These differences are well explained [here](https://www.davidegrayson.com/windev/msys2/) if you are interested.

## Get the tools

You can download MSYS2 [here](https://www.msys2.org/).
Once installed you open the shell from the start menu. This is a Linux-like terminal. Run the following line several times until is says no more updates:

`pacman -Syuu`

The close and re-open the shell and run this:

`pacman -S base-devel cmake gcc zlib`

## An example: aragorn

Aragorn is a program to detect tRNA genes and tmRNA genes in nucleotide sequences. It's a good example because there is no windows binary I could find and it's very easy to compile. You can download the source tarball [here](http://130.235.244.92/ARAGORN/Downloads/) and gunzip it. Then just follow the instructions on the aragorn page and run gcc like this:

`gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.38.c`

That's it. This gives you an executable aragorn.exe that you can now run. When you use it outside of MSYS2 or on another computer say, you will need to distribute it with the msys-2.0.dll file in the same folder.

Other c/c++ programs use a makefile for compiling so you will have to follow the specific instructions. Also if additional libraries are required for compilation you can install them with pacman. After compiling you can check which libraries your binary is linked to using the ldd command:

```
$ ldd aragorn.exe
        ntdll.dll => /c/WINDOWS/SYSTEM32/ntdll.dll (0x7ff816880000)
        KERNEL32.DLL => /c/WINDOWS/System32/KERNEL32.DLL (0x7ff816060000)
        KERNELBASE.dll => /c/WINDOWS/System32/KERNELBASE.dll (0x7ff813810000)
        msys-2.0.dll => /e/win_binaries/aragorn1.2.38/msys-2.0.dll (0x180040000)
```

## Links

* [MSYS2](https://www.msys2.org/)
* [Installing GCC & MSYS2](https://github.com/orlp/dev-on-windows/wiki/Installing-GCC--&-MSYS2)
* [MSYS2 explained better than I can](https://www.davidegrayson.com/windev/msys2/)
