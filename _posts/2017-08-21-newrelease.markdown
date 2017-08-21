---
layout: post
title:  "Version 0.8.1 released"
date:   2017-08-21 11:30:00
categories: dataexplore
tags: [releases,filtering]
---

Some fixes in version 0.8.1 include the ability to resample time series data with a datetime index, menu item for applying an element-wise function and a fix to a bug that made large tables draw slower.

## Improved filtering dialog

The filtering tool has now been changed so that you can add filters using widgets as well as (or instead of) using the string query. The example below shows multiple filters being chained together using the titanic dataset. There is no limit to the number of filters that can be used. If a string query is also entered it will be applied as well and chained to the first filter using the AND/OR logic selected. See [the wiki](https://github.com/dmnfarrell/pandastable/wiki/Filtering)

<div style="width: 400px;">
<a href="/img/filtering_dialog_example.png"><img src="/img/filtering_dialog_example.png" width="500px"></a>
</div>

## Undo changes

You can now undo the last applied operation though this doesn't yet apply to individual cell changes. Also there is only one level of undo at the moment. As always do NOT use this program as the primary storage place for valuable data and always keep the raw data.

## Links

Version 0.8.1 [release page](https://github.com/dmnfarrell/pandastable/releases/tag/v0.8.1).
