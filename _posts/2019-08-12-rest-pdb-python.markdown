---
layout: post
title:  "Accessing data from the PDB with Python"
date:   2019-08-12 12:40:00
categories: bioinformatics
tags: [pdb,python]
thumbnail: /img/rcsb-pdb.png
---

## Background

In the early days of the internet when programmers wanted to retrieve data from web pages automatically they would write code to read the page html directly and extract what they want. This is called 'screen scraping' and was a horribly inefficient and unreliable way to fetch data. If the web page was redesigned your code could break. If information was moved somewhere else you might have to write new code. The modern solution are so-called **RESTful** interfaces to web services. This just means you can write a piece of code to request specific data at predictable end points using an agreed set of rules. The rules are fixed and that data won't change so you always know what you are getting. You can retrieve the results as text, json or other formats. Many biological databases have REST services so that you can request data in bulk automatically.

The Protein Data Bank (PDB) supports RESTful services and is a good example to use here. REST services provide urls like normal web pages and can be viewed in a browser. So the following link can be accessed programmatically or in a browser: https://www.rcsb.org/pdb/rest/idStatus?structureId=4HHB. It returns xml or text format that you can parse. In the example below we want to retrieve  descriptions for a list of PDB IDs. This is much easier than going to each web page to check if you have hundreds.

## Code

This code requires pandas to run. It takes a list and converts it to a comma separated string. Then the request is formed and the request made using the `requests` module. This returns text which we finally parse into a pandas dataframe using `read_csv`. Websites will always have documentation about the form the rest string should take depending on the data you are looking for. In this case we have a `customReportColumns` field that we use to tell it the columns we want back.

```python
def get_pdb_descriptions(ids):
    """Takes a list of PDB codes and returns a dataframe with details for each one"""
    import pandas as pd
    import requests
    from io import StringIO

    pdbstr = ','.join(ids)
    #split the string up for readability
    pdbrest = 'http://www.rcsb.org/pdb/rest/customReport.xml?'\
    'pdbids=%s&customReportColumns=chainLength,uniprotRecommendedName,geneName,taxonomy'\
    '&service=wsfile&format=csv' %pdbstr
    #make the request
    r = requests.get(pdbrest)
    data = r.text
    #parse our returned csv data into a dataframe
    df = pd.read_csv(StringIO(data))    
    return df
  ```

We can then call the function as follows:

```python
ids = ['4esq','2h34']
result = get_pdb_descriptions(ids)
print (result)
```

The output is a dataframe as below. You can use the above method and the required URLs to get any csv type REST data back into a table with Python.

```
structureId chainId  chainLength         uniprotRecommendedName                geneName                    taxonomy
4ESQ       A          194  Serine/threonine-protein kinase PknH  pknH#Rv1266c#MTCY50.16  Mycobacterium tuberculosis
2H34       A          309  Serine/threonine-protein kinase PknE   pknE#Rv1743#MTCY28.05  Mycobacterium tuberculosis
2H34       B          309  Serine/threonine-protein kinase PknE   pknE#Rv1743#MTCY28.05  Mycobacterium tuberculosis
```

## Links

* [RCSB PDB RESTful Web Service interface](https://www.rcsb.org/pdb/software/rest.do)
* [UniProt REST access](https://www.uniprot.org/help/programmatic_access)
