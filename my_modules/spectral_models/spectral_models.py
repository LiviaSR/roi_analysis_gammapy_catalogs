#!/usr/bin/env python
# coding: utf-8

# [gammapy.modeling.models.SpectralModel](https://docs.gammapy.org/1.0/api/gammapy.modeling.models.SpectralModel.html#gammapy.modeling.models.SpectralModel)

# In[ ]:


import os
import sys
import importlib
path_my_modules = 'my_modules'
module_path = os.path.abspath(f'{path_my_modules}/config')
if module_path not in sys.path:
    sys.path.append(module_path)

import cfg
importlib.reload(cfg)


# In[ ]:


module_path = os.path.abspath(f'{path_my_modules}/{cfg.dir_utilities}')
if module_path not in sys.path:
    sys.path.append(module_path)

import utilities
importlib.reload(utilities)
from utilities import (
    mkdir_sub_directory, 
    write_tables_fits, 
    write_tables_csv, 
    load_catalogs_from_gammapy, 
    name_to_txt,
)


# In[2]:


from astropy import units as u
from gammapy.modeling.models import PowerLawSpectralModel, SkyModel

def sky_model_pl(
    index = 2,
    amplitude = 1e-12 * u.Unit("TeV-1 cm-2 s-1"),
    reference = 1 * u.Unit("TeV"),
    datasets_names = None
):
    """
    Returns a sky model with spectral model type: PowerLawSpectralModel
    >>>sky_model = sky_model_pl(
    index = 2,
    amplitude = 1e-12 * u.Unit("TeV-1 cm-2 s-1"),
    reference = 1 * u.Unit("TeV"),
    name = "pl",
)
    >>>print(sky_model)
    SkyModel

      Name                      : lp
      Datasets names            : None
      Spectral model type       : LogParabolaSpectralModel
      Spatial  model type       : 
      Temporal model type       : 
      Parameters:
        amplitude                     :   1.40e-12   +/- 6.7e-14 1 / (cm2 s TeV)
        reference             (frozen):      1.000       TeV         
        alpha                         :      1.577   +/-    0.03             
        beta                          :      0.233   +/-    0.01  
    """
    
    spectral_model = PowerLawSpectralModel(
        index=index, 
        amplitude=amplitude, 
        reference=reference,
    )
            
    if not datasets_names:
        sky_model = SkyModel(
            spectral_model = spectral_model,
        )
    else:
        sky_model = SkyModel(
            spectral_model = spectral_model, 
            name = f"{name_to_txt(datasets_names)}_{spectral_model.tag[1]}",
            datasets_names = datasets_names         
        )
    return sky_model


# In[1]:


from astropy import units as u
from gammapy.modeling.models import ExpCutoffPowerLawSpectralModel, SkyModel

def sky_model_ecpl(
    amplitude = 1e-12 * u.Unit("TeV-1 cm-2 s-1"),
    index = 2,
    lambda_= 0.1 * u.Unit("TeV-1"),
    reference = 10 * u.Unit("TeV"),
    alpha = 1.0,
    datasets_names = None
):
    """
    Returns a sky model with spectral model type: ExpCutoffPowerLawSpectralModel
    see: https://docs.gammapy.org/1.1/user-guide/model-gallery/spectral/plot_exp_cutoff_powerlaw.html
    
    >>>sky_model = sky_model_ecpl(
    amplitude = 1e-12 * u.Unit("TeV-1 cm-2 s-1"),
    index = 2,
    lambda_= 0.1 * u.Unit("TeV-1"),
    reference = 1 * u.Unit("TeV"),
    alpha = 1.0,
    name = "ecpl",
)
    >>>print(sky_model)
    SkyModel

      Name                      : lp
      Datasets names            : None
      Spectral model type       : LogParabolaSpectralModel
      Spatial  model type       : 
      Temporal model type       : 
      Parameters:
        amplitude                     :   1.40e-12   +/- 6.7e-14 1 / (cm2 s TeV)
        reference             (frozen):      1.000       TeV         
        alpha                         :      1.577   +/-    0.03             
        beta                          :      0.233   +/-    0.01  
    """
    
    spectral_model = ExpCutoffPowerLawSpectralModel(
        amplitude = amplitude,
        index = index,
        lambda_= lambda_,
        reference = reference,
        alpha = alpha,
    )
    if not datasets_names:
        sky_model = SkyModel(
            spectral_model = spectral_model,
        )    
    else:
        name = f"{name_to_txt(datasets_names)}_{spectral_model.tag[1]}"
        sky_model = SkyModel(
            spectral_model = spectral_model, 
            name = name,
            datasets_names = datasets_names
        )
    return sky_model


# In[ ]:





# In[ ]:


from astropy import units as u
from gammapy.modeling.models import LogParabolaSpectralModel, SkyModel

def sky_model_lp(
    alpha = 2.3,
    amplitude = 1e-12 * u.Unit("TeV-1 cm-2 s-1"),
    reference = 1 * u.Unit("TeV"),
    beta = 0.5,
    datasets_names = None
):
    """
    Returns a sky model with spectral model type: LogParabolaSpectralModel
    """
    
    spectral_model = LogParabolaSpectralModel(
        alpha = alpha,
        amplitude = amplitude,
        reference = reference,
        beta = beta,
    )
    if not datasets_names:
        sky_model = SkyModel(
            spectral_model = spectral_model,
        )    
    else:
        name = f"{name_to_txt(datasets_names)}_{spectral_model.tag[1]}"
        sky_model = SkyModel(
            spectral_model=spectral_model, 
            name= name,
            datasets_names = datasets_names
        )
    return sky_model


# In[ ]:


from astropy import units as u
import matplotlib.pyplot as plt
from gammapy.modeling.models import BrokenPowerLawSpectralModel, Models, SkyModel

def sky_model_bpl(
    index1=1.5,
    index2=2.5,
    amplitude=1e-12 * u.Unit("TeV-1 cm-2 s-1"),
    ebreak=1 * u.Unit("TeV"),
    datasets_names = None
):
    """
    Returns a sky model with spectral model type: BrokenPowerLawSpectralModel
    """
    
    spectral_model =BrokenPowerLawSpectralModel(
        index1=index1,
        index2=index2,
        amplitude=amplitude,
        ebreak=ebreak,
    )
    if not datasets_names:
        sky_model = SkyModel(
            spectral_model = spectral_model,
        )    
    else:
        name = f"{name_to_txt(datasets_names)}_{spectral_model.tag[1]}"
        sky_model = SkyModel(
            spectral_model=spectral_model, 
            name= name,
            datasets_names = datasets_names 
        )
    return sky_model


# In[ ]:


import pandas as pd 
from gammapy.datasets import FluxPointsDataset
from gammapy.catalog import CATALOG_REGISTRY
import os
import sys

# format_csv = ".csv"
format_fits = ".fits"

sed_type="e2dnde"

lst=[]
ds_lst = []

def get_source_data(dict_lhaaso_tevc = None, catalog_tags = None, path_dir = None):
    
    i_range = range(len(dict_lhaaso_tevc.keys()))
    for i in i_range:

        LHAASO_name = list(dict_lhaaso_tevc.keys())[i]
        LHAASO_id = LHAASO_name.replace(" ", "")
        
        j_range = range(len(dict_lhaaso_tevc[LHAASO_name]))
        for j in j_range:

            catalog_src = []

            pf_on = []
            src_on = []

            ds_j=[]

            source_name=dict_lhaaso_tevc[LHAASO_name][j]
            source_id = source_name.replace(" ", "")
            
            k_range = range(len(catalog_tags))
            for k in k_range: 
                
                catalog_tag = catalog_tags[k]
                catalog=CATALOG_REGISTRY.get_cls(catalog_tag)()
                
                try:
                    
                    src=catalog[source_name]
                    src_on.append(src.data)
                    catalog_src.append(catalog_tag)
                    
                    ds = FluxPointsDataset(
                        data=src.flux_points, 
                        name=catalog_tag
                    )
                    
                    ds_j.append(ds)
                    pf_on.append(catalog_tag)

                    table = ds.data.to_table(
                        sed_type=sed_type, 
                        formatted=True
                    )

                    file_name = f'{LHAASO_id}_{source_id}_{catalog_tag}{format_fits}'
                    path_os = os.path.abspath(
                        os.path.join(
                            f"{path_dir}/{file_name}"
                        )
                    )
                    
                    
                    if path_os not in sys.path:
                        sys.path.append(path_os)

                    #table.write(f"{path_os}{format_csv}",format='ascii.ecsv', overwrite=True)
                    table.write(f"{path_os}",format='fits', overwrite=True)
                    
                except:
                    pass

                lst_k=[LHAASO_name, source_name, catalog_src, pf_on, ds_j, src_on]
            lst.append(lst_k)
            ds_lst.append(ds_j)

    df = pd.DataFrame(lst, columns =['LHAASO', 'TeV Conterpart', 'Catalog', 'Flux Points', 'ds', 'src']) 
    df.to_csv(f"{path_dir}/flux_points.csv", index = True )
    return df, ds_lst

