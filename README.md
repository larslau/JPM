# Joint Propagation Model
The Joint propagation model (JPM) is a model developed for deriving a harmonized quantification of tau PET data across different tracers. The JPM simultaneously use head-to-head data and anchor point data to estimate equations for mapping regional SUVRs to a harmonized scale named *CenTauRs*. 

The JPM is implemented in the R and C++ using the [Template Model Builder](https://kaskr.github.io/adcomp/_book/Introduction.html) framework.

> [!NOTE]
> If you are using this code, please reference
> 
> Leuzy, A., Raket, L.L., et al. "Roadmap for harmonizing tau PET in Alzheimer's disease: the Joint Propagation model". *Alzheimer's & Dementia* (2024).

Example
--------------------
Below is an illustration of how JPM can produce harmonized CenTauR units from SUVR obtained from different tracers. Simulations can be run using [example.R](example.R) which also illustrates how to use the JPM code.

![](img/CTR_mapping.png)
![](img/CTR_eval.png)
