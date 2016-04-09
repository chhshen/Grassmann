----------------------
Citation
----------------------

Please cite the following article when using this source code:

  M. Harandi, C. Sanderson, C. Shen and B. Lovell. ``Dictionary Learning and Sparse Coding on Grassmann Manifolds: An Extrinsic Solution’’.
  Int. Conference on Computer Vision (ICCV), 2013.



```
  @InProceedings{Harandi2013ICCV,
    author = {Harandi, Mehrtash and Sanderson, Conrad and Shen, Chunhua and Lovell, Brian C.},
    title = {Dictionary Learning and Sparse Coding on Grassmann Manifolds: An Extrinsic Solution},
    booktitle = {The IEEE International Conference on Computer Vision (ICCV)},
    month = {December},
    pages={3120-3127},
    doi={10.1109/ICCV.2013.387},
    ISSN={1550-5499},
    year = {2013}
  }
```

Or,

```
@article{Harandi2015IJCV,
   author    = "M. Harandi and  R. Hartley and  C. Shen and  B. Lovell and  C. Sanderson",
   title     = "Extrinsic methods for coding and dictionary learning on {G}rassmann manifolds",
   journal   = "International Journal of Computer Vision",
   volume    = "114",
   number    = "2",
   year      = "2015",
 }
```



- DOI: 10.1109/ICCV.2013.387

You can obtain a copy of this article via:
http://www.cv-foundation.org/openaccess/content_iccv_2013/papers/Harandi_Dictionary_Learning_and_2013_ICCV_paper.pdf



----------------------
License
----------------------

The source code is provided without any warranty of fitness for any purpose.
You can redistribute it and/or modify it under the terms of the
GNU General Public License (GPL) as published by the Free Software Foundation,
either version 3 of the License or (at your option) any later version.
A copy of the GPL license is provided in the "GPL.txt" file.



----------------------
Instructions and Notes
----------------------

This code uses some functions from the SPAMS (SPArse Modeling Software) toolbox (http://spams-devel.gforge.inria.fr/). For the ease of the user, we have incorporated the required functions from the SPAMS into the package.
As such, there is hopefully no need to install the SPAMS toolbox. If you face any error regarding SPAMS functions, please download and install the SPAMS package.

To execute the code, run the `Run_gSC.m` for sparse coding (labelled dictionary is required) and `Run_gSC_dic.m` for dictionary learning and sparse coding. We also provide the Cambridge Hand Gesture dataset (Kim etal. CVPR'07) for evaluation purposes.
The evaluation procedure is to use the first 80 samples of each class as probe data and the remaining samples as gallery set. To cite the Cambridge hand gesture dataset use the following item

```
  @InProceedings{Kim_2007_CVPR,
				 title={Tensor canonical correlation analysis for action classification},
				 author={Kim, Tae-Kyun and Wong, Kwan-Yee Kenneth and Cipolla, Roberto},
				 booktitle={IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
				 pages={1--8},
				 year={2007},
				 organization={IEEE}
 }
```




