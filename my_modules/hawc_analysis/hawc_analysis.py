#!/usr/bin/env python
# coding: utf-8

# In[ ]:





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


from gammapy.modeling.models import Models
from gammapy.datasets import Datasets
from gammapy.modeling.models import PowerLawSpectralModel
def get_dataset(source_name):
    # sky model HAWC J1825-134
    # https://arxiv.org/pdf/2012.15275.pdf
    if source_name == "HAWC J1825-134":
        sky_model = spec.sky_model_pl(
            index=2.28,
            amplitude="4.2e-15 TeV-1 cm-2 s-1",
            reference=18 * u.TeV,
        datasets_names = source_name
        )

    # sky model HAWC J1825-138
    # https://arxiv.org/pdf/2012.15275.pdf
    if source_name == "HAWC J1825-138":
        sky_model = spec.sky_model_ecpl(
            amplitude=2.7e-14 * u.Unit("cm-2 s-1 TeV-1"),
            index=2.02,
            lambda_=1/27 * u.Unit("TeV-1"),
            reference=18 * u.TeV,
        datasets_names = source_name
        )

    # sky model HAWC J1826-128
    # https://arxiv.org/pdf/2012.15275.pdf
    if source_name == "HAWC J1826-128":
        sky_model = spec.sky_model_ecpl(
            amplitude=2.7e-14 * u.Unit("cm-2 s-1 TeV-1"),
            index=1.2,
            lambda_=(1/24) * u.Unit("TeV-1"),
            reference=18 * u.TeV,
        datasets_names = source_name
        )

    # sky model eHWC J1825-134
    # https://arxiv.org/pdf/1909.08609.pdf
    if source_name == "eHWC J1825-134":
        sky_model = spec.sky_model_ecpl(
            amplitude=2.12e-13 * u.Unit("cm-2 s-1 TeV-1"),
            index=2.12,
            lambda_= (1/61) * u.Unit("TeV-1"),
            reference=10 * u.TeV,
        datasets_names = source_name
        )

    # sky model eHWC J1907+063
    # https://arxiv.org/pdf/1909.08609.pdf
    if source_name == "eHWC J1907+063":
        sky_model = spec.sky_model_lp(
             alpha=2.46,
            amplitude="0.95e-13 cm-2 s-1 TeV-1",
            reference=10 * u.TeV,
            beta=0.11,
        datasets_names = source_name
        )

    # sky model eHWC J2019+368
    # https://arxiv.org/pdf/1909.08609.pdf
    if source_name == "eHWC J2019+368":
        sky_model = spec.sky_model_lp(
            alpha=2.08,
            amplitude="0.45e-13 cm-2 s-1 TeV-1",
            reference=10 * u.TeV,
            beta=0.26,
        datasets_names = source_name
        )
            # print()
    if source_name == "2HWC J1825-134":
        sky_model = spec.sky_model_pl(
            index=2.58,
            amplitude="138.0e-15 cm-2 s-1 TeV-1",
            reference=7 * u.TeV,
            datasets_names = source_name
        )
        # print() 
            
    table = table_to_SED_format(cfg.path_fp_HAWC, utl.name_to_txt(source_name))  
    return utl.ds_fp_from_table_fp(table = table, sky_model = sky_model, source_name=source_name)


# In[ ]:





# In[ ]:





# In[ ]:





# In[2]:


from astropy.coordinates import SkyCoord

def get_dict():
    '''
    Dictionary of the HAWC sources (dict keys: Source name based on J2000 coordinates) information: "position": Right ascension (in degrees), Declination (in degrees) 
    '''
    unit_deg = cfg.unit_deg
    return {
    "2HWC J1825-134": {
        "position":  SkyCoord(276.46,-13.40, unit=unit_deg), # 2HWC catalogue       
    }, 
    "HAWC J1825-138": {
        "position":  SkyCoord(276.38,-13.86, unit=unit_deg), # https://arxiv.org/pdf/2012.15275.pdf
    },
    "HAWC J1826-128": {
        "position":  SkyCoord(276.50,-12.86, unit=unit_deg) # https://arxiv.org/pdf/2012.15275.pdf
    },
    "HAWC J1825-134": {
        "position":  SkyCoord(276.44,-13.42, unit=unit_deg) # https://arxiv.org/pdf/2012.15275.pdf
    },
    "eHWC J1825-134": {
        "position":  SkyCoord(276.40,-13.37, unit=unit_deg) # https://arxiv.org/pdf/1909.08609.pdf
    },
    "eHWC J1907+063": {
        "position":  SkyCoord(286.91,6.32, unit=unit_deg) # https://arxiv.org/pdf/1909.08609.pdf
    }, 
    "eHWC J2019+368": {
        "position":  SkyCoord(304.95,36.78, unit=unit_deg) # https://arxiv.org/pdf/1909.08609.pdf
    }, 
}


# In[1]:


from astropy.table import Table
from astropy import units as u
from pathlib import Path
def table_to_SED_format(path, file_name):
    '''
    Normalization Representation
    The SED format is a flexible specification for representing one-dimensional spectra 
    (distributions of amplitude vs. energy).
    
    '''
    toBool = {'True':True,'False':False}

    file_path = Path(f'{path}/{file_name}{cfg.format_csv}') 
    table = Table.read(file_path,format='ascii', delimiter=',', comment='#')
    
    
#     display(table)

    table['is_ul'] = True

    for index, row in enumerate(table['col7']):
        if table['col7'][index] == 'True':
            table['is_ul'][index] = True
        else:
            table['is_ul'][index] = False
            
    table.rename_column('col1', 'e_ref')
    table['e_ref'].unit = u.TeV
  
    table.rename_column('col2', 'e_min')
    table['e_min'].unit = u.TeV

    table.rename_column('col3', 'e_max')
    table['e_max'].unit = u.TeV
    
    table.rename_column('col4', 'e2dnde')
    table['e2dnde'].unit = u.Unit("TeV cm-2 s-1")
        
    table['e2dnde_ul'] = table['e2dnde']
#     table['e2dnde_ul'].unit = u.Unit("TeV cm-2 s-1")
    
    table.rename_column('col5', 'e2dnde_errp')
    table['e2dnde_errp'].unit = u.Unit("TeV cm-2 s-1")
    table.rename_column('col6', 'e2dnde_errn')
    table['e2dnde_errn'].unit = u.Unit("TeV cm-2 s-1")
    # tablesky_model['is_ul'].dtype('bool')
    table.meta["SED_TYPE"] = "e2dnde"
    
    l1 = table['e_min']
    l2 = table['e_max']

    if(set(l1) == set(l2)):
        table_out = table['e_ref', 'e2dnde', 'e2dnde_errp', 'e2dnde_errn', 'e2dnde_ul', 'is_ul'].copy()
    else:
        table_out = table['e_ref', 'e_min', 'e_max', 'e2dnde', 'e2dnde_errp', 'e2dnde_errn', 'e2dnde_ul', 'is_ul'].copy()

    display(table_out)
    
    return table_out


# In[1]:


path_fp_HAWC = "/home/born-again/Documents/GitHub/CTA_projects/flux_points_outside_gammapy_catalogs/HAWC"
from gammapy.modeling.models import Models
from gammapy.datasets import Datasets

def get_HAWC_tables_datasets(dict_HAWC):
    models = Models()
    for source_index, source_name in enumerate(list(dict_HAWC.keys())):
         # sky model HAWC J1825-134
        # https://arxiv.org/pdf/2012.15275.pdf
        if source_name == "HAWC J1825-134":
            sky_model = sky_model_pl(
                index=2.28,
                amplitude="4.2e-15 TeV-1 cm-2 s-1",
                reference=18 * u.TeV,
                datasets_names = source_name
            )
            # print(sky_model)
            models.append(sky_model.copy(name = sky_model.name,  datasets_names = source_name))
            # print(models)

        # sky model HAWC J1825-138
        # https://arxiv.org/pdf/2012.15275.pdf
        if source_name == "HAWC J1825-138":
            sky_model = sky_model_ecpl(
                amplitude=2.7e-14 * u.Unit("cm-2 s-1 TeV-1"),
                index=2.02,
                lambda_=1/27 * u.Unit("TeV-1"),
                reference=18 * u.TeV,
                datasets_names = source_name
            )
            # print(sky_model)
            models.append(sky_model.copy(name = sky_model.name,  datasets_names = source_name))
            # print(models)

        # sky model HAWC J1826-128
        # https://arxiv.org/pdf/2012.15275.pdf
        if source_name == "HAWC J1826-128":
            sky_model_HAWCJ1826_128 = sky_model_ecpl(
                amplitude=2.7e-14 * u.Unit("cm-2 s-1 TeV-1"),
                index=1.2,
                lambda_=(1/24) * u.Unit("TeV-1"),
                reference=18 * u.TeV,
                datasets_names = source_name
            )
            # print(sky_model_HAWCJ1826_128)
            models.append(sky_model_HAWCJ1826_128.copy(name = sky_model_HAWCJ1826_128.name,  datasets_names = source_name))
            # print(models)

        # sky model eHWC J1825-134
        # https://arxiv.org/pdf/1909.08609.pdf
        if source_name == "eHWC J1825-134":
            sky_model = sky_model_ecpl(
                amplitude=2.12e-13 * u.Unit("cm-2 s-1 TeV-1"),
                index=2.12,
                lambda_= (1/61) * u.Unit("TeV-1"),
                reference=10 * u.TeV,
                datasets_names = source_name
            )
            # print(sky_model)
            models.append(sky_model.copy(name = sky_model.name,  datasets_names = source_name))
            # print(models)

        # sky model eHWC J1907+063
        # https://arxiv.org/pdf/1909.08609.pdf
        if source_name == "eHWC J1907+063":
            sky_model = sky_model_lp(
                 alpha=2.46,
                amplitude="0.95e-13 cm-2 s-1 TeV-1",
                reference=10 * u.TeV,
                beta=0.11,
                datasets_names = source_name
            )
            # print(sky_model)
            models.append(sky_model.copy(name = sky_model.name,  datasets_names = source_name))
            # print(models)

        # sky model eHWC J2019+368
        # https://arxiv.org/pdf/1909.08609.pdf
        if source_name == "eHWC J2019+368":
            sky_model = sky_model_lp(
                alpha=2.08,
                amplitude="0.45e-13 cm-2 s-1 TeV-1",
                reference=10 * u.TeV,
                beta=0.26,
                datasets_names = source_name
            )
            # print(sky_model)
            models.append(sky_model.copy(name = sky_model.name,  datasets_names = source_name))
            # print(models)

            # print()
        if source_name == "2HWC J1825-134":
            sky_model = sky_model_pl(
                index=2.58,
                amplitude="138.0e-15 cm-2 s-1 TeV-1",
                reference=7 * u.TeV,
                datasets_names = source_name
            )
            # print() 
            models.append(sky_model.copy(name = sky_model.name,  datasets_names = source_name))
            # print(models)

    tables = []
    datasets = Datasets()
    for source_index, source_name in enumerate(list(dict_HAWC.keys())):
        print(source_index, source_name)
        table = table_to_SED_format(path_fp_HAWC, utl.utl.name_to_txt(source_name))
        tables.append(table)
        dataset = utilities.ds_fp_from_table_fp(table = table,  sky_model= models[source_index], source_name=source_name)
        datasets.append(dataset)    
    return tables, datasets, models


# In[ ]:





# In[ ]:




