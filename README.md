## Safety-critical traffic control using CLB
A reproduction of the result in the [[paper]](https://arxiv.org/pdf/2301.04833.pdf):  
Zhao, Chenguang, Huan Yu, and Tamas G. Molnar. "Safety-critical traffic control by connected automated vehicles." Transportation research part C: emerging technologies 154 (2023): 104230  

### Results of sudden brake of the leading HDV
* nominal controller <br><img src="simulation/nominal.gif" alt="nominal.gif" width="70%" />  
* STC controller <br><img src="simulation/STC.gif" alt="STC.gif" width="70%" />  
### Run the code
language: Julia v1.6.7
```shell
Julia
using IJulia
notebook()
```
Then play code in the demo.ipynb

