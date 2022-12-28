---
layout: post
title:  "Scrape dynamic tables in Python with Playwright"
date:   2022-12-29 11:13:00
categories: python
tags: [python]
thumbnail: /img/web_scraping_dalle.png
---

## Background

<div style="width: 350px; float:right;">
 <a href="/img/web_scraping_dalle.png"> <img src="/img/web_scraping_dalle.png" width="300px"></a>
</div>

Sometimes it's necessary to scrape a website or some pages that contain elements generated dynamically, often via javascript. This means you can't always get all the html content in the page directly from the url. You may have to simulate interaction in the browser like clicking a button on a form and then waiting to get the content. This can be done with the **Playwright** library which was created to support automated testing of websites. Note that some websites won't like you scraping their pages in some cases and you should be careful about making many requests in a short time from a website.

## Example

In this example we want to extract the elements in a table that is generated when a form in the page is filled in. The resulting table is paginated and has buttons at the bottom to go to the next page of results. We want to extract all these results automatically and save the html table elements in a list. We can then parse these table rows with BeautifulSoup. The way it's done here is to launch a browser from playwright and actually load the page, then wait a bit for all the bits to load. The example site here is fictional. We want to search for items between certain dates. These are filled in by finding the form element names on the web page first. Use the console for this.

In this case the inputs are called `input#DateFrom` and `input#DateFrom`. We then loop over the number of results pages (if we know this) and get the table element by name (with id `page_table`) using `page.inner_html`. This html is appended to a list. Finally the 'Next' button is clicked to go to the next page in the table before we start the loop again. After this we close the browser.

```python
from playwright.sync_api import Playwright, sync_playwright
from bs4 import BeautifulSoup

def get_search_results(start="17/12/22",end="18/12/22"):
    """Fetch a dynamic html table from a form"""
    res = []
    html = 'test'
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True,slow_mo=50)
        page = browser.new_page()
        page.goto('https://examplewebsite/com/search')
        # Fill in the dates
        page.fill("input#DateFrom", start)
        page.fill("input#DateTo", end)

        for i in range(5):
            #wait to load table fully
            page.wait_for_timeout(1000)
            # Extract the HTML of the results
            html = page.inner_html("#page_table")
            res.append(html)
            #get next page in table
            page.locator('"Next"').click()

        # Close the browser
        browser.close()
    return res
```

This gives us a list of lists for each page with the table row elements in them as below:

```
<td class="th_name" style="border-left: 1px solid rgb(211, 211, 211);"><div class="item_table" onclick="location.href='/show.php?x/y/z"><span>surname'"><div class="p_name" style="text-align:left;"><a class="showdn-link" href="/show.php?x/y/z"><span>surname</span>, Paul</a></div><div class="clear"></div></div></td>
```

We can then just extract the data we want from the html with bs4. In this case we are trying to get the urls from the table. They are in the `<a>` element with the class `showdn-link`.

```python
rows = get_search_results()
links = []
for html in rows:
    soup = BeautifulSoup(html, 'html.parser')
    rows = soup.find_all('a',class_="showdn-link")
    for a in rows:    
        links.append(pref+a['href'])
```

This is just one way to use Playwright without even having to 'crawl' the site since we are looking for a specific table result.

## Links

* [Playwright](https://playwright.dev/python/docs/intro)
* [Beautiful Soup Documentation](https://www.crummy.com/software/BeautifulSoup/bs4/doc/)
