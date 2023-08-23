<p align="center">

  <h1 align="center">
    SVIndelGenotyper
  </h1>

  <p align="center">
   <a href="https://github.com/stjude/SVIndelGenotyper" target="_blank">
     <img alt="Status"
          src="https://img.shields.io/badge/status-active-success.svg" />
   </a>
   <a href="https://github.com/stjude/SVIndelGenotyper/issues" target="_blank">
     <img alt="Github Issues"
          src="https://img.shields.io/github/issues/stjude/SVIndelGenotyper"  />
   </a>
   <a href="https://github.com/stjude/SVIndelGenotyper/pulls"  target="_blank">
     <img alt="Pull Requests"
          src="https://img.shields.io/github/issues-pr/stjude/SVIndelGenotyper"  />
   </a>
   <a href="https://github.com/stjude/SVIndelGenotyper/blob/main/LICENSE" target="_blank">
     <img alt="License: Apache-2.0"
          src="https://img.shields.io/github/license/saltstack/salt" />
   </a>
  </p>


  <p align="center">
   SVIndelGenotyper was employed for analyzing the error profiles of indels and structural variants in deep sequencing data.    
   <br />
   <a href="#"><strong>Explore the docs »</strong></a> 
   <br />
   <a href="#"><strong>Read the paper »</strong></a>
   <br /> 
   <br />
   <a href="https://github.com/stjude/SVIndelGenotyper/issues/new?assignees=&labels=&template=feature_request.md&title=Descriptive%20Title&labels=enhancement">Request Feature</a>
    | 
   <a href="https://github.com/stjude/SVIndelGenotyper/issues/new?assignees=&labels=&template=bug_report.md&title=Descriptive%20Title&labels=bug">Report Bug</a>
   <br />
    ⭐ Consider starring the repo! ⭐
   <br />
  </p>
</p>

---
## Quick Start


### Installation
SVIndelGenotyper depends on [python>=3.9.9](https://www.python.org/downloads/)

```
> git clone https://github.com/stjude/SVIndelGenotyper.git
> cd SVIndelGenotyper
> ./make_install.sh
```

## Usage 
SVIndelGenotyper has 2 subcommands:
* ```indel``` genotype indels for sequencing data (WGS/WES/RNASeq)
* ```sv``` genotype SVs for sequencing data (WGS)

```
> module load python/3.9.9
> python SVIndelGenotyper.py indel testIndel.tsv
> python SVIndelGenotyper.py sv testSV.tsv
```

### Contact
* Pandurang.Kolekar[AT]stjude.org
* Xiaotu.Ma[AT]stjude.org

---
#### COPYRIGHT 
Copyright 2023 St. Jude Children's Research Hospital

#### LICENSE
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
