---
layout: post
title:  "Convert a multi-sample VCF to a pandas DataFrame"
date:   2020-11-28 12:20:00
categories: bioinformatics
tags: [python,genomics,ngs]
thumbnail: /img/multi_sample_vcf_igv.png
---

## Background

Here is some code I wrote to convert a vcf file with many samples into a table format. This was done to make the calls for easy to read at a given site. Reading a multi sample vcf is tortuous. The vcf is read in using pyVCF and for each record (a site with calls) the calls for each sample are parsed with the depth values (in the calldata object). We create a list for each and at the end convert all into a dataframe with the right column names. This code does assume your vcf has the FORMAT fields included as follows. It doesn't generalise to other formats. Also this code was only tested on a vcf made with `bcftools call` on a bacterial genome alignment.

```
LT708304.1	1057	.	A	G	999	.	DP=5494;ADF=1,3028;ADR=0,1869;AD=1,4897;VDB=0.206626;SGB=24.3331;RPB=1;MQB=1;MQSB=1;BQB=1;MQ0F=0;AC=46;AN=46;DP4=1,0,3034,1873;MQ=60	GT:PL:DP:SP:ADF:ADR:AD	1:255,0:78:0:0,44:0,33:0,77	1:255,0:114:0:1,67:0,46:1,113	1:255,0:96:0:0,57:0,39:0,96	1:255,0:98:0:0,59:0,39:0,98	1:255,0:64:0:0,34:0,30:0,64	1:255,0:116:0:0,75:0,41:0,116	1:255,0:114:0:0,71:0,43:0,114	1:255,0:56:0:0,30:0,25:0,55	1:255,0:111:0:0,58:0,53:0,111	1:255,0:131:0:0,83:0,47:0,130	1:255,0:236:0:0,139:0,96:0,235	1:255,0:146:0:0,98:0,48:0,146	1:255,0:115:0:0,74:0,41:0,115	1:255,0:84:0:0,56:0,28:0,84	1:255,0:141:0:0,94:0,47:0,141
```

## Code

```python
import pandas as pd
import vcf
from gzip import open as gzopen

def vcf_to_dataframe(vcf_file):
    """
    Convert a multi sample vcf to dataframe. Records each samples FORMAT fields.
    Args:
        vcf_file: input multi sample vcf
    Returns: pandas DataFrame
    """

    ext = os.path.splitext(vcf_file)[1]
    if ext == '.gz':
        file = gzopen(vcf_file, "rt")
    else:
        file = open(vcf_file)
    vcf_reader = vcf.Reader(file,'r')
    res=[]
    cols = ['sample','REF','ALT','mut','DP','ADF','ADR','AD','chrom','var_type','sub_type','start','end','QUAL']
    i=0
    for rec in vcf_reader:
        #if i>50:
        #    break
        x = [rec.CHROM, rec.var_type, rec.var_subtype, rec.start, rec.end, rec.QUAL]
        for sample in rec.samples:
            if sample.gt_bases == None:
              #no call
                mut=''
                row = [sample.sample, rec.REF, sample.gt_bases, mut, 0,0,0,0]
            elif rec.REF != sample.gt_bases:
                mut = str(rec.end)+rec.REF+'>'+sample.gt_bases
                cdata = sample.data
                row = [sample.sample, rec.REF, sample.gt_bases, mut, cdata[2], cdata[4] ,cdata[5], cdata[6]] + x
            else:
                #call is REF
                mut = str(rec.end)+rec.REF              
                cdata = sample.data
                row = [sample.sample, rec.REF, sample.gt_bases, mut, cdata[2], cdata[4] ,cdata[5], cdata[6]] + x

            res.append(row)
    res = pd.DataFrame(res,columns=cols)
    res = res[~res.start.isnull()]
    return res

df = vcf_to_dataframe('test.vcf.gz')  
```

Now instead of having to read the vcf file and look at each site as above, we can look at the calls for all samples at this site in a table format as below. This is the same information as the vcf excerpt given above for position 1057:

|     sample | REF | ALT |     mut |  DP |      ADF |     ADR |       AD |      chrom | var_type | sub_type | start  | end    | QUAL  |
|-----------:|----:|----:|--------:|----:|---------:|--------:|---------:|-----------:|----------|----------|--------|--------|-------|
|   13-11594 |   A |   G | 1057A>G |  78 |  [0, 44] | [0, 33] | [0, 77]  | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
|  14-MBovis |   A |   G | 1057A>G | 114 |  [1, 67] | [0, 46] | [1, 113] | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
|   15-11643 |   A |   G | 1057A>G |  96 |  [0, 57] | [0, 39] | [0, 96]  | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
|   17-11662 |   A |   G | 1057A>G |  98 |  [0, 59] | [0, 39] | [0, 98]  | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
|  17-MBovis |   A |   G | 1057A>G |  64 |  [0, 34] | [0, 30] | [0, 64]  | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
| 182-MBovis |   A |   G | 1057A>G | 116 |  [0, 75] | [0, 41] | [0, 116] | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
| 19-11957   | A   | G   | 1057A>G | 114 | [0, 71]  | [0, 43] | [0, 114] | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
| 19-MBovis  | A   | G   | 1057A>G | 56  | [0, 30]  | [0, 25] | [0, 55]  | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
| 22-12200   | A   | G   | 1057A>G | 111 | [0, 58]  | [0, 53] | [0, 111] | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
| 23-MBovis  | A   | G   | 1057A>G | 131 | [0, 83]  | [0, 47] | [0, 130] | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
| 24-MBovis  | A   | G   | 1057A>G | 236 | [0, 139] | [0, 96] | [0, 235] | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
| 25-MBovis  | A   | G   | 1057A>G | 146 | [0, 98]  | [0, 48] | [0, 146] | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
| 26-12883   | A   | G   | 1057A>G | 115 | [0, 74]  | [0, 41] | [0, 115] | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
| 26-MBovis  | A   | G   | 1057A>G | 84  | [0, 56]  | [0, 28] | [0, 84]  | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |
| 27-MBovis  | A   | G   | 1057A>G | 141 | [0, 94]  | [0, 47] | [0, 141] | LT708304.1 | snp      | ts       | 1056.0 | 1057.0 | 999.0 |


## Links

* [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format)
* [pyVCF](https://pyvcf.readthedocs.io/en/latest/index.html)
* [snpgenie](https://github.com/dmnfarrell/snpgenie)
