```
«Покоряя вершину за вершиной, мы открываем полные великолепия и интереса просторы,
но мы не видим нашей цели, мы не видим горизонта;
вдали возвышаются еще более высокие вершины,
которые откроют тем, кто взойдет на них,
еще более широкие перспективы и углубят чувство,
истинность которого подчеркивается каждым достижением науки, что
“Велики дела Господни”».

Инаугурационная речь на ежегодном собрании Британской  Научной Ассоциации проф. 
сэра Джозефа Джона Томсона
26 августа 1909 г.
```

# dynastatflex
Python library for modeling dynamic (now only static) of flexible mechanisms. Useful for robotics application.

There is pipline for future development and equations derivation: [MODELING_WITH_DAMPING.md](logsec/MODELING_WITH_DAMPING.md). To correct view formulas you should use logseq.

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
- - - -










