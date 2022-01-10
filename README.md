# VCFStatistics
A method for visualizing VCF data on a circos plot (data is generated in a format easy to use with the Circos software)

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#requirements">Requirements</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>

<!-- requirements -->
## Requirements

These scripts have been tested with Python 2.7 and 3 and should work with either.
These scripts requires a VCF file (with a depth field for certain metrics).  The VCF file can be compressed with gzip or bgzip.

<!-- usage -->
## Usage

1) Find metrics for each individual in each window

      python VCF.stats.windowBased.v1.0.py -vcf file.vcf -win 10000 -depth 3 > statsfile.out

      To see the usage and get further information: python VCF.stats.windowBased.v1.0.py -h

2) Summarize data for groups and possibly output in a format that can be used with Circos software to visualize.

      python PartitionData.LsalSpecific.py -file statsfile.out -ind ind1,ind3 -avg no -out depth -win 10000 -chr no > output.txt

      To see the usage and get further information: python PartitionData.LsalSpecific.py -h
      
Optionally, identify regions that differ by missing genotypes between different groups (used for a specific project)

      python SexSpecific.LsalSpecific.py -file statsfile.out -ind1 ind1,ind3 -ind2 ind2,ind4 -win 10000 -fold 2
      
      To see the usage and get further information: python SexSpecific.LsalSpecific.py -h

<!-- license -->
## License 

Distributed under the MIT License.
