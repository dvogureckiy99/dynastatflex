# dynastatflex
Python library for modeling dynamic (now only static) of flexible mechanisms. Useful for robotics application.

## Dependencies
If you want use MATLAB SPACAR package to compare the results you need MATLAB. 
- - - -
Was tested with 2023a version.

## Installation
First need to setup conda environment:
```
conda env create --file dynastatflex.yml -n dynastatflex
```
If you using VSC, then reload it. Choose dynastatflex environment.
- - - - 
If you want to update conda env with .yml file type in Command Prompt:
```
conda env update --file dynastatflex.yml --prune
```
## Examples
![alt text](images/momentcenterbeam_SPACARmatch.png)
- - - -
![alt text](images/momentcenterbeam_SPACARmatch2.png)

## Referencies:
[1] Huber G, Wollherr D and Buss M (2021) A Concise and Geometrically Exact Planar Beam Model for Arbitrarily Large Elastic Deformation Dynamics. Front. Robot. AI 7:609478. doi: 10.3389/frobt.2020.609478. link=[https://www.frontiersin.org/articles/10.3389/frobt.2020.609478/full](https://www.frontiersin.org/articles/10.3389/frobt.2020.609478/full)
