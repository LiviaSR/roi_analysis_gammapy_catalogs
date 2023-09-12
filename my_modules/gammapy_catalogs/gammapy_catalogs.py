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


module_path = os.path.abspath(f'{path_my_modules}/{cfg.dir_hawc_analysis}')
if module_path not in sys.path:
    print(sys.path.append(module_path))
    sys.path.append(module_path)

import hawc_analysis as hawc
importlib.reload(hawc)


# In[ ]:


module_path = os.path.abspath(f'{path_my_modules}/{cfg.dir_lhaaso_analysis}')
if module_path not in sys.path:
    print(sys.path.append(module_path))
    sys.path.append(module_path)

import lhaaso_analysis as lhaaso
importlib.reload(lhaaso)


# In[ ]:


module_path = os.path.abspath(f'{path_my_modules}/{cfg.dir_hess_analysis}')
if module_path not in sys.path:
    print(sys.path.append(module_path))
    sys.path.append(module_path)

import hess_analysis as hess
importlib.reload(hess)


# In[4]:


def get_datasets_flux_points_outside_gammapy(region_of_interest):
    '''
    Select a catalog subset (only sources within a region of interest)
    '''
    path_datasets = utl.get_path_datasets(region_of_interest)
    
    datasets_counterparts = Datasets() # global datasets object
    models_counterparts = Models()  # global models object
    souces_names = []
    n_counterparts = 0 # number of counterparts
    n_flux_points = 0 # number of flux points tables
    
    e_ref_min = region_of_interest["e_ref_min"] 
    e_ref_max = region_of_interest["e_ref_max"]
    
    roi_pos = region_of_interest["position"]
    radius_roi = region_of_interest["radius_roi"]

    dict_HAWC = hawc.get_dict()
    dict_LHAASO = lhaaso.get_dict()
    for index, source_name in enumerate(list(dict_HAWC.keys())):
        source_pos = dict_HAWC[source_name]["position"]
        if roi_pos.separation(source_pos) <= radius_roi:
            ds = hawc.get_dataset(source_name)
            if any([e_ref_min !=  None, e_ref_max !=  None]):
                    ds = cut_flux_points_energy_range(ds, e_ref_min, e_ref_max)
            souces_names.append(source_name)        
            counterpart_model = ds.models[0]
            models_counterparts.append(counterpart_model)
            datasets_counterparts.append(ds)
            
            n_flux_points+=1
            n_counterparts+=1 
            
    for index, source_name in enumerate(list(dict_LHAASO.keys())):
        source_pos = dict_LHAASO[source_name]["position"]
        if roi_pos.separation(source_pos) <= radius_roi:
            ds = lhaaso.get_dataset(source_name)
            if any([e_ref_min !=  None, e_ref_max !=  None]):
                    ds = cut_flux_points_energy_range(ds, e_ref_min, e_ref_max)
            souces_names.append(source_name) 
            counterpart_model = ds.models[0]
            models_counterparts.append(counterpart_model)
            datasets_counterparts.append(ds)

            n_flux_points+=1
            n_counterparts+=1
            
    if datasets_counterparts:
        datasets_counterparts.models = models_counterparts
        # To save datasets and models
        utl.write_datasets_models(datasets_counterparts,region_of_interest, "counterparts_outside_gammapy")
            
        print(f"Total number of counterparts: {n_counterparts}")
        print(f"Total number of flux points tables: {n_flux_points}")
    
        return souces_names, datasets_counterparts, models_counterparts
    else: print("No counterparts in the ROI")


# In[ ]:


import numpy as np
from gammapy.estimators import FluxPoints

def cut_flux_points_energy_range(dataset, e_ref_min, e_ref_max):
    
    if all([e_ref_min ==  None, e_ref_max ==  None]):
            raise Exception(f"Sorry, there is a error: e_ref_min is None and e_ref_max is None)") 
    
    flux_points = dataset.data
    models = dataset.models[0]      
    ds_name = dataset.name
    
    if e_ref_min != None:
        mask_energy = np.zeros(len(flux_points.to_table()), dtype=bool)

        for m, e_ref in enumerate(flux_points.energy_ref):
            if e_ref >= e_ref_min:
                mask_energy[m] = True

        flux_points_mask = flux_points.to_table()[mask_energy]
        flux_points = FluxPoints.from_table(flux_points_mask)
    
    if e_ref_max != None:
        mask_energy = np.zeros(len(flux_points.to_table()), dtype=bool)

        for m, e_ref in enumerate(flux_points.energy_ref):
            if e_ref <= e_ref_max:
                mask_energy[m] = True

        flux_points_mask = flux_points.to_table()[mask_energy]
        flux_points = FluxPoints.from_table(flux_points_mask)     
        
    return FluxPointsDataset(models = models, data = flux_points, name = ds_name)


# In[ ]:


# from gammapy.datasets import FluxPointsDataset
from astropy.coordinates import SkyCoord
from astropy import units as u

import pickle

def create_catalogs_region_of_interest(region_of_interest):
    """
    Gets catalogs subset (only sources within the radius of the region of interest)
    """
    
    source_catalogs = utl.load_catalogs_from_gammapy()
    
    source_position = region_of_interest["position"] 
    radius_roi = region_of_interest["radius_roi"] 
        
    catalogs_roi = []
    catalogs_no_counterparts = []
    numbers_catalogs_roi = 0
    
    for catalog in source_catalogs:        
        # Selects only sources within the region of interest. 
        mask_roi = source_position.separation(catalog.positions) < radius_roi 
        
        if len(catalog[mask_roi].table):
            catalogs_roi.append(catalog[mask_roi])
            numbers_catalogs_roi += 1
        else:
            catalogs_no_counterparts.append(f"{catalog.tag}: {catalog.description}")
    
            
    if numbers_catalogs_roi:
        utl.pickling_catalog_roi(catalogs_roi, region_of_interest)
        print(f"\n{numbers_catalogs_roi} catalogs with sources within the region of interest:", end = "\n\n")
        for catalog in catalogs_roi:
            print(f"{catalog.tag}: {catalog.description}")
            display(catalog.table)
    else:
        print("No catalogs with sources in the region of interest!", end = "\n\n")

    if numbers_catalogs_roi and len(catalogs_no_counterparts):
        print("Catalogs without sources within the region of interest:", end = "\n\n")
        for index, catalog_no_counterpart in enumerate(catalogs_no_counterparts):                            
            print(catalog_no_counterpart)

    return catalogs_roi


# In[ ]:


from gammapy.datasets import FluxPointsDataset
from astropy.coordinates import SkyCoord
from astropy import units as u
from gammapy.modeling.models import SkyModel, Models
from gammapy.datasets import Datasets

def get_datasets_flux_points_gammapy(region_of_interest):
    '''
    Select a catalog subset (only sources within a region of interest)
    '''
    
    path_datasets = utl.get_path_datasets(region_of_interest)
    try:
        catalogs_roi = utl.unpickling_catalog_roi(region_of_interest)
    except:
        catalogs_roi = create_catalogs_region_of_interest(region_of_interest)
        
    datasets_counterparts = Datasets() # global datasets object
    models_counterparts = Models()  # global models object
    counterparts = [] # global sources object
    
    n_counterparts = 0 # number of counterparts
    n_flux_points = 0 # number of flux points tables
    
    e_ref_min = region_of_interest["e_ref_min"] 
    e_ref_max = region_of_interest["e_ref_max"]
    
    #############################
    if region_of_interest["name"] == 'LHAASO J1839-0545':
        for cat_tag in ["gamma-cat", "hgps"]:
            hess.get_dataset_HESSJ1837069(counterparts, datasets_counterparts, models_counterparts, cat_tag)
    #############################
        
    for catalog in catalogs_roi:
        cat_tag = catalog.tag
        for counterpart in catalog:
            n_counterparts+=1   
            counterpart_name = counterpart.name            
            try:
                flux_points = counterpart.flux_points

                counterpart_spectral_model = counterpart.spectral_model()
                spectral_model_tag = counterpart_spectral_model.tag[0]
                spectral_model_tag_short = counterpart_spectral_model.tag[1]
                
                if cat_tag == 'gamma-cat' or cat_tag == 'hgps':
                    ds_name = f'{counterpart_name}: {cat_tag}'
                else: ds_name = counterpart_name
                    
                file_name = utl.name_to_txt(ds_name)
                
                counterpart_model = SkyModel(
                    name = f"{file_name}_{counterpart_spectral_model.tag[1]}",
                    spectral_model = counterpart_spectral_model,
                    datasets_names=ds_name
                )
        
                ds = FluxPointsDataset(
                    models = counterpart_model,
                    data = flux_points, 
                    name =  ds_name   
                )
                
                if any([e_ref_min !=  None, e_ref_max !=  None]):
                    ds = cut_flux_points_energy_range(ds, e_ref_min, e_ref_max)
                
                n_flux_points+=1
                models_counterparts.append(counterpart_model)  # Add the counterpart_model to models()
                
                counterparts.append(counterpart)
                datasets_counterparts.append(ds)
                
            except Exception as error:
                # By this way we can know about the type of error occurring
                print(f'The error is: ({counterpart_name}) {error}') 
            
    datasets_counterparts.models = models_counterparts
    # To save datasets and models
    utl.write_datasets_models(datasets_counterparts,region_of_interest, "counterparts_gammapy")
            
    print(f"Total number of counterparts: {n_counterparts}")
    print(f"Total number of flux points tables: {n_flux_points}")
    return counterparts, datasets_counterparts, models_counterparts

