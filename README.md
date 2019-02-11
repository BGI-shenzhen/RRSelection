

# RRSelection
<b>RRSelection:  A new simple and efficient software to detect selection region analysis based Variant Call Format</b>
  
	RRSelection: A Llinkage disequilibrium method to detect selection region across population VCF

###  1) Install
------------

<b> [Download](https://github.com/BGI-shenzhen/RRSelection/archive/v0.85.tar.gz) </b>

</br>
Method1 For <b>linux/Unix</b> and <b> macOS </b>
<pre>
        git clone https://github.com/BGI-shenzhen/RRSelection.git
	cd RRSelection ;chmod 755 configure; ./configure;
        make;
        mv RRSelection  bin/;    #     [rm *.o]
</pre>

**Note:** If fail to link,try to <b>re-install</b> the libraries [**_zlib_**](https://zlib.net/)

Method2 For <b>linux/Unix</b> and <b> macOS </b>
<pre>
        tar -zxvf  RRSelectionXXX.tar.gz
        cd RRSelectionXXX;
        cd src;
        make ; make clean                            # or [sh make.sh]
        ../bin/RRSelection
</pre>
**Note:** If fail to link,try to <b>re-install</b> the libraries [**_zlib_**](https://zlib.net/)

###  2) Example
------------

see more detailed Usage in the <b>[Documentation](https://github.com/BGI-shenzhen/RRSelection/blob/master/Manual.txt)</b>

* 1) Calculate sliding windows mean RR for one or two population,and give out the selection region. also give out the whole genome  RR plot figure. 
<pre>
      # 1)  For all samples in one population 
               ./bin/RRSelection   -InVCF SNP.vcf.gz  -OutPut OutPrefix
      # 2)  For same samples in one population
	       ./bin/RRSelection   -InVCF SNP.vcf.gz  -OutPut OutPrefix  -SubGroup  subgroup.list  # subgroup.list is the sample name of this population
      # 3)  For Tow  population
               ./bin/RRSelection   -InVCF SNP.vcf.gz  -OutPut OutPrefix  -SubGroup  subgroup.list  #  PopID : sample name list
</pre>

* 2) see the result  [OutPrefix.winRR.gz OutPrefix.selection.gz] and [OutPrefix.png OutPrefix.pdf]. ALso Run the perl script to regain the beautiful picture

<pre>
         perl     PlotRRSele.pl    -inFile   OutPrefix.winRR.gz  -output OutPrefix
</pre>
 
    

###  3) Introduction
------------
To detect the selection region is the most important and most common analysis in the population resequencing. Here we introduce a new software :<b>RRSelection</b>, a simple-efficient software to detect the selection region analysis based Variant Call Format. Sliding whole genome windows to calculate every region mean R^2 for one or two population,  and pick out the high-chained region (one population) or the region with the greatest difference (two populations), which is regarded the selection region according to the top  distribution of measure MeanRR(Z-test Pvalue).


* <b> Parameter description</b>
```php
        Usage: RRSelection  -InVCF  <in.vcf.gz>  -OutPut <outPrefix>

               -InVCF      <str>     Input SNP VCF Format
               -OutPut     <str>     OutPut sliding stat mean r^2 Result

               -SubGroup   <str>     one/two sub-group Sample List File,-h for more help
               -Windows    <int>     Sliding windows bin (kb),MaxDis between two pairwise SNP[300]
               -Step       <float>   Step ratio(0,1] of windows,1:NoOverlap [0.2]
               -Masked     <int>     Masked windows when the SNP Number too low[10]

               -MAF        <float>   Min minor allele frequency filter [0.05]
               -Het        <float>   Max ratio of het allele filter [0.88]
               -Miss       <float>   Max ratio of miss allele filter [0.25]

               -Pvalue     <float>   T-test Pvalue to pick out selection region[0.005]
               -KeepR                Keep Rscript used to modify and plots

	       -help                 See more help [hewm2008 Beta v0.85]

```


###  4) Results
------------
The following  is the format of the result output file header .  and the Figure is no showed here.
<pre>
#Chr    Start   End     Mean_r^2_cul    Sum_r^2_cul     Count_cul       Mean_r^2_wild   Sum_r^2_wild    Count_wild      MeanRRDiff(cul-wild)    ZScore  Pvalue
##Group[cul], MeanRR:0.245096   SD:0.0529981    Effective windows Count:30
##Group[wild], MeanRR:0.247118  SD:0.0631814    Effective windows Count:30
##Diff MeanRR[cul-wild], Mean:-0.00202159       SD:0.0276476    Effective windows Count:30
Tu2     0       300000  0.1779  68833.0295      386921  0.1034  60456.6791      584831  0.0745  2.77    0.002822
Tu2     60000   360000  0.0782  12577.9534      160803  0.0734  21013.7645      286430  0.0049  0.25    0.401158
Tu2     120000  420000  0.1223  26006.7348      212590  0.0877  31059.4803      354331  0.0347  1.33    0.092056
...
</pre>
###  5) Discussing
------------
- [:email:](https://github.com/BGI-shenzhen/RRSelection) hewm2008@gmail.com / hewm2008@qq.com
- join the<b><i> QQ Group : 125293663</b></i>

######################swimming in the sky and flying in the sea #############################
