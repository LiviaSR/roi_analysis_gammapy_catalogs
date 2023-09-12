#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import sys
import importlib
path_my_modules = "/home/born-again/Documents/GitHub/CTA_projects/my_modules"
module_path = os.path.abspath(f'{path_my_modules}/config')
if module_path not in sys.path:
    sys.path.append(module_path)

import cfg
importlib.reload(cfg)


# In[1]:


import os
import sys
import importlib

module_path = os.path.abspath(f'{path_my_modules}/{cfg.dir_spectral_models}')
if module_path not in sys.path:
    sys.path.append(module_path)

import spectral_models as spec
importlib.reload(spec)

module_path = os.path.abspath(f'{path_my_modules}/{cfg.dir_utilities}')
if module_path not in sys.path:
    sys.path.append(module_path)

import utilities as utl
importlib.reload(utl)


# In[ ]:


from gammapy.catalog import CATALOG_REGISTRY 
from gammapy.modeling.models import SkyModel

def get_dataset_HESSJ1837069(counterparts, datasets_counterparts, models_counterparts, cat_tag = "gamma-cat"):
    source_name = "HESS J1837-069"

    catalog = CATALOG_REGISTRY.get_cls(cat_tag)()
    source = catalog[source_name]
    if cat_tag == "gamma-cat":
        spectral_model = source.sky_model().spectral_model
    else: 
        spectral_model = spec.sky_model_pl().spectral_model
    
    datasets_names = f'{source_name}: {cat_tag}'
    name=f"{utl.name_to_txt(datasets_names)}_{spectral_model.tag[1]}"
    table = source.flux_points.to_table(sed_type=cfg.sed_type_e2dnde)
    sky_model = SkyModel(
        spectral_model=spectral_model,   
        name=name,
        datasets_names = datasets_names
    )
    dataset = utl.ds_fp_from_table_fp(table = table, sky_model = sky_model, source_name=datasets_names)
    counterparts.append(source) 
    datasets_counterparts.append(dataset)
    models_counterparts.append(sky_model)

