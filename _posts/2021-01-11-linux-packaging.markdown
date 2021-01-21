---
layout: post
title:  "Linux application packaging and universal formats"
date:   2021-01-11 12:43:00
categories: software
tags: [linux,snapcraft]
thumbnail: img/linux-packaging.png
---

## Linux package management

<div style="width: 250px; float:right;">
 <a href="/img/linux-packaging.png"> <img src="/img/linux-packaging.png" width="200px"></a>
</div>

Linux is an amazing open source operating system used universally across industry, in science, on embedded hardware and on servers. The core system and kernel are typically packaged up into 'distributions' for end users. So we have Ubuntu, Fedora, Debian and endless others. These are really variants on the same thing but some are specialised for desktop or server use (or something more specific like security). Actually there are far too many distributions but only a few are widely used. They add a lot of extra functionality for users and modern desktop linux systems are in many ways now better than their commercial counterparts.

For application developers there is a problem with all this. Distros use different packaging systems for their applications (e.g. dnf, apt, eopkg). These often (always?) rely on linking to shared libraries known to be present on that specific version of the OS. If you want your app to work it must be available in the package ecosystem and getting it there can be onerous. Firstly it can be complicated to make a package from your application and then you have to maintain it for each version of the distro. Or you can get it put into the official repositories and maintained by someone else but that means a lot more hurdles. This can be enormously off putting and not worth the effort for niche or little known applications. For popular tools this won't be a problem so much as it will be curated by the Ubuntu, Fedora or whatever developers package applications.

Linus Torvalds, the creator of Linux, has aptly summarised the problem at a conference in 2014:

<blockquote>
We basically don't make binaries for Linux. Why? Because binaries for Linux desktop applications is a major f*ing pain in the ass. Right. You don't make binaries for Linux. You make binaries for Fedora 19, Fedora 20, maybe there's even like RHEL 5 from ten years ago, you make binaries for debian stable. - Linus Torvalds
</blockquote>

See more quotes in the explanation of this issue on the [AppImage github page](https://github.com/AppImage/AppImageKit).

## A solution

Ideally the developer would build their own package which will then run on basically any Linux distribution regardless of the distribution or version. This would contain with it most of the dependencies needed to run. It would be built in such a way that it's likely to run on current and recent stable versions of the distros. This way you can get the package straight from the developer with fewer hurdles. Solutions like this are now available in several forms.

## Snap, Flatpak and AppImage

The three main such universal formats are Snaps, Flatpak and AppImage. The general idea is that the app runs in a container with limited access to the host system (a 'sandbox'). AppImage is the oldest and most straightforward solution in some ways. All the dependencies are bundled together using a build recipe and the file (.AppImage) is shared as a single binary which is simply downloaded and run. They have limited security features so could be potentially risky to use. Though the files can be signed. The packages work on a large number of systems. There is a website, [AppImageHub](https://www.appimagehub.com/) for hosting them though you can download them from anywhere.

The Snap format is made by canonical (makers of Ubuntu) and the files are called [snaps](https://en.wikipedia.org/wiki/Snap_(package_manager)). The snap file format is a single compressed filesystem using the SquashFS format with the extension `.snap`. Snaps are braoder in scope and work for packaging desktop apps, server tools, IoT apps and even system services. Snaps only work on distributions that use the init system. You can install snaps from the [snap store](https://snapcraft.io/store).

Flatpak is also container formats which run on a wide number of distributions. There is a central repository at [Flathub](flathub.org). It is only meant for desktop apps. Flatpak and snaps need a service running to use the packages, AppImages don't. Snaps have a far superior security model than the others.

How easy is it to put your apps into these formats? Well it depends. **Snapcraft** is a tool for developers to package their snaps. I only have experience using Snapcraft and AppImage (a bit). Snapcraft is easy to get started with and now has decent documentation and a forum where the snapcraft developers provide help. Early adoption was frustrating because the tool changed so much and so rapidly. Snapcraft builds on virtual machines using multipass and can use a bit of disk space in the process. If your application is complex you may have to do a of work to get it working. AppImage building seems to have a less steep learning curve.

## Adoption

I don't think any of these formats are yet very widely adopted compared to the traditional package managers - yet. AppImage seems to intended to allow small developers especially package up their apps and bypass the distro gatekeepers. Snap is the most ambitious, with the weight of Canonical behind it. They have met resistance in the Linux community in adopting them or gaining traction outside of the Ubuntu ecosystem. But the format has such a wide scope, like on servers and devices, that even if it fails on the desktop it will not be going away. I think any of these formats are good thing for developers who can't (or don't want to) get their apps into the official repositories. They let you give users the latest version of your software directly in a sandboxed format isolated from the host system. All of these formats are improving over time and if you are a developer looking to package your tool, give at least one of them a try. These are likely to be the future of application package delivery, eventually..

## Links

* [snapcraft home](https://snapcraft.io/build)
* [snapcraft forum](https://forum.snapcraft.io/)
* [flatpak](https://flatpak.org/)
* [AppImage](https://github.com/AppImage/AppImageKit)
* [Comparison: Snap vs Flatpak vs AppImage](https://linuxhint.com/snap_vs_flatpak_vs_appimage/)
* [Snap vs. Flatpak vs. AppImage: Know The Differences, Which is Better](https://www.fosslinux.com/42410/snap-vs-flatpak-vs-appimage-know-the-differences-which-is-better.htm)
