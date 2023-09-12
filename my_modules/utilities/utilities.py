#!/usr/bin/env python
# coding: utf-8

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


module_path = os.path.abspath(f'{path_my_modules}/{cfg.dir_hawc_analysis}')
if module_path not in sys.path:
    print(sys.path.append(module_path))
    sys.path.append(module_path)

import hawc_analysis as hawc


# In[ ]:


module_path = os.path.abspath(f'{path_my_modules}/{cfg.dir_cta_simulation}')
if module_path not in sys.path:
    print(sys.path.append(module_path))
    sys.path.append(module_path)

import cta_simulation as cta


# In[ ]:


from gammapy.modeling import Fit
from gammapy.datasets import Datasets
from gammapy.modeling.models import Models

def fit_Datasets(datasets, sky_model):
    model_name = sky_model.name
    datasets = Datasets(datasets)
    sky_model = sky_model.copy(name = model_name)
        
    datasets.models = sky_model
    # print(datasets)
    fitter = Fit()
    result_fit = fitter.run(datasets=datasets)
    print(result_fit.parameters.to_table())
    print(result_fit.total_stat)
    
#     models = Models(sky_model.copy(name= model_name)) 
#     file_path = utl.get_path_models(region_of_interest)
#     # To save only the models
#     models.write(f"{file_path}/{model_name}.yaml", model_name)
    region_of_interest = create_region_of_interest()

    write_datasets_models(datasets, region_of_interest, model_name)
    return sky_model


# In[ ]:


# // In your script file
def getVarFromFile(filename):
    import imp
    f = open(filename)
    global data
    data = imp.load_source('data', filename, f)
    f.close()
    return data


# In[ ]:


def write_datasets_models(datasets,region_of_interest, directory_name, path_datasets = None):
    
    if path_datasets == None:
        path_datasets = get_path_datasets(region_of_interest)        
    
    path_datasets, path_file = mkdir_sub_directory(str(path_datasets), directory_name)
    datasets.write(filename=f"{path_file}/datasets{cfg.format_yaml}", filename_models=f"{path_file}/models{cfg.format_yaml}", overwrite=True)
    return


# In[ ]:


def read_datasets_models(region_of_interest, directory_name, path_datasets = None):
    if path_datasets == None:
        path_datasets = get_path_datasets(region_of_interest)  
    path_datasets, path_file = mkdir_sub_directory(str(path_datasets), directory_name)
    return Datasets.read(filename=f"{path_file}/datasets{cfg.format_yaml}", filename_models=f"{path_file}/models{cfg.format_yaml}")


# In[ ]:


def create_file_name(pulsar_info, irf_name, skymodel_name): 
    """
    """
    return f"{name_to_txt(pulsar_info['name'])}_irf_{irf_name}_{skymodel_name}"


# In[ ]:


from gammapy.catalog import CATALOG_REGISTRY 

def load_catalogs_from_gammapy():
    """
    Load all available source catalogs in gammapy.catalog package into a list
    
    >>> load_catalogs_from_gammapy()
    
    Source catalogs in Gammapy: 8

    (catalog index: 0) SourceCatalogGammaCat:
        name: gamma-cat
        description: An open catalog of gamma-ray sources
        sources: 162

    (catalog index: 1) SourceCatalogHGPS:
        name: hgps
        description: H.E.S.S. Galactic plane survey (HGPS) source catalog
        sources: 78

    (catalog index: 2) SourceCatalog2HWC:
        name: 2hwc
        description: 2HWC catalog from the HAWC observatory
        sources: 40

    (catalog index: 3) SourceCatalog3FGL:
        name: 3fgl
        description: LAT 4-year point source catalog
        sources: 3034

    (catalog index: 4) SourceCatalog4FGL:
        name: 4fgl
        description: LAT 8-year point source catalog
        sources: 6659

    (catalog index: 5) SourceCatalog2FHL:
        name: 2fhl
        description: LAT second high-energy source catalog
        sources: 360

    (catalog index: 6) SourceCatalog3FHL:
        name: 3fhl
        description: LAT third high-energy source catalog
        sources: 1556

    (catalog index: 7) SourceCatalog3HWC:
        name: 3hwc
        description: 3HWC catalog from the HAWC observatory
        sources: 65
    """
    source_catalogs = []
    print (f"Source catalogs in Gammapy: {len(CATALOG_REGISTRY)}\n")
    for index, catalog in enumerate(CATALOG_REGISTRY):
        #  FITS files are loaded
        catalog_cls = CATALOG_REGISTRY.get_cls(catalog.tag)()
        source_catalogs.append(catalog_cls)
        print(f"(catalog index: {index}) {source_catalogs[index]}")

    return source_catalogs


# In[ ]:


from astropy.coordinates import SkyCoord

def set_source_info():
    """
    Sets the source info into a dictionary

    Parameters
    ----------
    source_name : str
        source name based on J2000 coordinates
        
    source_RA : float  
        right ascension (in degrees) of the source position
        
    source_dec : float
        declination (in degrees) of the source position
         
    Returns
    -------
    source_info : dict 
        dictionary with the source info (name and position)
    
    Example
    -------
    >>> source_name = "LHAASO J1825-1326"  
    >>> source_RA = 276.45* u.Unit("deg")  
    >>> source_dec = -13.45* u.Unit("deg") 
    >>> set_source_info(source_name, source_RA, source_dec)
    {'name': 'LHAASO J1825-1326',
    'position': <SkyCoord (ICRS): (ra, dec) in deg
    (276.45, -13.45)>}
    """
    data = getVarFromFile("set_analysis.dat")
    source_name = data.source_name
    pos_ra = data.pos_ra*u.Unit(cfg.unit_deg)
    pos_dec = data.pos_dec*u.Unit(cfg.unit_deg)

#     if any([source_RA.unit !=  cfg.unit_deg, source_dec.unit !=  cfg.unit_deg]):
#         raise Exception("Sorry, there is a error: celestial coordinates (RA, dec.) units is not in degrees") 
    
    return  {
        'name': source_name,
        'position': SkyCoord(pos_ra, pos_dec) 
    }


# In[ ]:


from astropy import units as u

def create_region_of_interest():
    """
    Creates the region of interest
    
    Parameters
    ----------
    source_info : dict 
        dictionary with the source info (name and position)
        
    radius_roi : float
        the maximum angle (in degrees) of separation between the source and its counterpart
         
    Returns
    -------
    region_of_interest : dict 
        dictionary with the region of interest info (source name, source position, angle of separation and roi name)
    
    Example
    -------
    >>> source_name = "LHAASO J1825-1326"  
    >>> source_RA = 276.45  
    >>> source_dec = -13.45 
    >>> source_info = set_source_info(source_name, source_RA, source_dec)
    >>> radius_roi = 1
    >>> create_region_of_interest(source_info, radius_roi)
    
    {'name': 'LHAASO J1825-1326',
     'position': <SkyCoord (ICRS): (ra, dec) in deg
     (276.45, -13.45)>,
     'radius_roi': <Quantity 1. deg>}
     
    """
    source_info = set_source_info()

    data = getVarFromFile("set_analysis.dat")
    radius_roi = data.radius_roi*u.Unit(cfg.unit_deg)
    e_ref_min = data.e_ref_min
    e_ref_max = data.e_ref_max
    
    if radius_roi.unit !=  cfg.unit_deg:
        raise Exception("Sorry, there is a error: radius_roi unit is not in degrees") 
        
    if e_ref_min is not None and e_ref_max is not None:
        if e_ref_min.unit != e_ref_max.unit:
            raise Exception(f"Sorry, there is a error: units is not iquals ({e_ref_min.unit} != {e_ref_max.unit})") 
            
        if e_ref_max.value <= e_ref_min.value:
            raise Exception(f"There is a error: e_ref_max ({e_ref_max}) <= e_ref_min ({e_ref_min})") 
            
    region_of_interest = source_info.copy()
    region_of_interest["radius_roi"] = radius_roi
    region_of_interest["e_ref_min"] = e_ref_min
    region_of_interest["e_ref_max"] = e_ref_max
    region_of_interest["roi_name"] = create_roi_name(source_info['name'], radius_roi, e_ref_min, e_ref_max)
    
#     print("**************************************************", end = "\n\n")
#     print(f"Region of interest:\n")
#     print(f'Source name: {region_of_interest['name']}')
#     print(f'Source position (ra, dec) in deg: {region_of_interest['position'].ra.deg, region_of_interest['position'].dec.deg}')
#     print(f"Radius in deg: {radius_roi.value}")
#     print(f"Energy ref min in deg: {radius_roi.value}")
    
#     print(f"Radius in deg: {radius_roi.value}", end = "\n\n**************************************************\n")
    
    return region_of_interest


# In[ ]:





# In[ ]:


def create_roi_name(source_name, radius_roi, e_ref_min = None, e_ref_max = None): 
    """
    """
    source_name = f"{name_to_txt(source_name)}"
    radius_name = f"_roi_{name_to_txt(str(radius_roi.value))}{name_to_txt(str(radius_roi.unit))}"

    if e_ref_min is not None:
        e_ref_min_name = f"_e_ref_min_{name_to_txt(str(e_ref_min.value))}{e_ref_min.unit}"
    else:
        e_ref_min_name = ""
        
    if e_ref_max is not None:
        e_ref_max_name = f"_e_ref_max_{name_to_txt(str(e_ref_max.value))}{e_ref_max.unit}"
    else:
        e_ref_max_name = ""   
    
    return(f"{source_name}{radius_name}{e_ref_min_name}{e_ref_max_name}") 


# In[ ]:


def create_data_frame_counterparts(region_of_interest):
    import pandas as pd 
    df_columns=[]
    df = pd.DataFrame()
    
    catalogs_roi = unpickling_catalog_roi(region_of_interest)

    for catalog in catalogs_roi:
        cat_tag = catalog.tag

        for counterpart in catalog:
            try: 
                flux_points = counterpart.flux_points
                flux_points_table = "Yes"
            except:
                flux_points_table = "No"
                                     
            sep = counterpart.position.separation(region_of_interest['position']).deg
            pos_dec = counterpart.position.dec.deg
            pos_ra = counterpart.position.ra.deg
            if cat_tag != 'gamma-cat' and cat_tag != 'hgps':
                name = counterpart.name
            else:
                name = f"{counterpart.name} ({cat_tag})"
            df_column = [name,"{:.2f}".format(pos_ra) , "{:.2f}".format(pos_dec), "{:.2f}".format(sep), flux_points_table]
            df_columns.append(df_column)

    df = pd.DataFrame(df_columns, columns = ['Source name', 'RA(deg)', 'dec.(deg)', 'Sep.(deg)', 'Flux points']) 
    df = df.reset_index(drop = True)
    df.index.name = 'Source index'
    return df


# In[ ]:


def get_dict_sources_gammapy(sources_gammapy, datasets_gammapy):
    dict_sources_gammapy = {}
    for index, (source, dataset) in enumerate(zip(sources_gammapy, datasets_gammapy)):
        dict_sources_gammapy[dataset.name] = {'position': source.position}
    return dict_sources_gammapy


# In[ ]:


from gammapy.modeling.models import Models
from gammapy.datasets import Datasets

def joint_datasets(datasets_gammapy, models_gammapy, datasets_outside_gammapy, models_outside_gammapy):
    datasets_roi = Datasets()
    models_roi = Models()
    for index, (dataset, model) in enumerate(zip(datasets_gammapy, models_gammapy)):
        datasets_roi.append(dataset)
        models_roi.append(model)
        print(index, dataset.name)
    for index_, (dataset, model) in enumerate(zip(datasets_outside_gammapy, models_outside_gammapy)):
        datasets_roi.append(dataset)
        models_roi.append(model)
        print(index_+index+1, dataset.name)
    datasets_roi.models = models_roi
    return datasets_roi, models_roi


# In[ ]:


import pandas as pd 
def get_dict_data_frame_roi(sources_gammapy, datasets_gammapy, datasets_roi, region_of_interest):
    
    df_columns=[]
    df = pd.DataFrame()
    dict_sources_gammapy = get_dict_sources_gammapy(sources_gammapy, datasets_gammapy)
    dict_HAWC = hawc.get_dict()
    dict_pulsars = cta.get_dict_pulsars()

    dict_roi = {}
    columns = ['Source name', 'RA(deg)', 'dec.(deg)', 'Sep.(deg)']
    for index, source_name in enumerate(datasets_roi.names):
        roi_pos = region_of_interest["position"]
        radius_roi = region_of_interest["radius_roi"]

        if source_name in list(dict_sources_gammapy.keys()):
            source_pos = dict_sources_gammapy[source_name]["position"]
            pos_dec = source_pos.dec.deg
            pos_ra = source_pos.ra.deg
            sep = source_pos.separation(roi_pos).deg
            df_column = [source_name,"{:.2f}".format(pos_ra) , "{:.2f}".format(pos_dec), "{:.2f}".format(sep)]
            df_columns.append(df_column)
            dict_roi[source_name] = {
                'position': source_pos,
                'separation':sep }

        if source_name in list(dict_HAWC.keys()):
            source_pos = dict_HAWC[source_name]["position"]
            pos_dec = source_pos.dec.deg
            pos_ra = source_pos.ra.deg
            sep = source_pos.separation(roi_pos).deg
            df_column = [source_name,"{:.2f}".format(pos_ra) , "{:.2f}".format(pos_dec), "{:.2f}".format(sep)]
            df_columns.append(df_column)
            dict_roi[source_name] = {
                'position': source_pos,
                'separation':sep }
    #     else: 
    #         source_pos = dict_LHAASO[name]["position"]
    #         sep = source_pos.separation(roi_pos).deg
    #         print(name, sep)

    for index, source_name in enumerate(list(dict_pulsars.keys())):
            source_pos = dict_pulsars[source_name]["position"]
            pos_dec = source_pos.dec.deg
            pos_ra = source_pos.ra.deg
            sep = source_pos.separation(roi_pos)
            if sep <= radius_roi:
                sep = sep.deg
                df_column = [source_name,"{:.2f}".format(pos_ra) , "{:.2f}".format(pos_dec), "{:.2f}".format(sep)]
                df_columns.append(df_column)
                dict_roi[source_name] = {
                    'position': source_pos,
                    'separation':sep }

    df = pd.DataFrame(df_columns, columns = columns) 
    df = df.reset_index(drop = True)
    return dict_roi, df


# In[ ]:


from gammapy.datasets import FluxPointsDataset
from astropy.coordinates import SkyCoord
from astropy import units as u
from gammapy.modeling.models import SkyModel, Models
from gammapy.datasets import Datasets

def get_flux_points_datasets(region_of_interest):
    '''
    Select a catalog subset (only sources within a region of interest)
    '''
    
    # Creates the directories to save the flux points tables 
#     directory_file = f'{cfg.dir_flux_points_tables}/{region_of_interest["roi_name"]}'
#     path_analysis, path_file = mkdir_sub_directory(cfg.dir_analysis, directory_file)
    path_file = get_path_tables(region_of_interest)
    catalogs_roi = unpickling_catalog_roi(region_of_interest)

    datasets_counterparts = Datasets() # global datasets object
    models_counterparts = Models()  # global models object
    counterparts = [] # global sources object
    
    n_counterparts = 0 # number of counterparts
    n_flux_points = 0 # number of flux points tables

    for catalog in catalogs_roi:
        cat_tag = catalog.tag
        for counterpart in catalog:
            n_counterparts+=1   
            counterpart_name = counterpart.name            
            try: 
                counterpart_flux_points = counterpart.flux_points
                n_flux_points+=1
                counterpart_spectral_model = counterpart.spectral_model()
                spectral_model_tag = counterpart_spectral_model.tag[0]
                spectral_model_tag_short = counterpart_spectral_model.tag[1]
        
                if cat_tag != 'gamma-cat' and cat_tag != 'hgps':
                    ds_name = f"{counterpart_name}"
                else:
                    ds_name = f"{counterpart_name}: {cat_tag}"
                     
                file_name = name_to_txt(ds_name)
    
                counterpart_model = SkyModel(
                    name = f"{file_name}_{counterpart_spectral_model.tag[1]}",
                    spectral_model = counterpart_spectral_model,
                    datasets_names=ds_name
                )
                models_counterparts.append(counterpart_model)  # Add the counterpart_model to models()
        
                ds = FluxPointsDataset(
                    models = counterpart_model,
                    data = counterpart_flux_points, 
                    name =  ds_name   
                )
                counterparts.append(counterpart)
                datasets_counterparts.append(ds)
                
                table = ds.data.to_table(sed_type = cfg.sed_type_e2dnde, formatted = True)

                # Writes the flux points table in the csv/fits format
                write_tables_csv(table, path_file, file_name)
                write_tables_fits(table, path_file, file_name)
                
            except Exception as error:
                # By this way we can know about the type of error occurring
                print(f"The error is: ({counterpart_name}) {error}") 
                            
    print(f"Total number of counterparts: {n_counterparts}")
    print(f"Total number of flux points tables: {n_flux_points}")
    return counterparts, datasets_counterparts, models_counterparts


# In[ ]:


# from gammapy.datasets import FluxPointsDataset
from astropy.coordinates import SkyCoord
from astropy import units as u

import pickle

def get_catalogs_region_of_interest(source_catalogs, region_of_interest):
    """
    Gets catalogs subset (only sources within the region of interest)
    """
    
    source_position = region_of_interest['position'] 
    radius_roi = region_of_interest["radius_roi"] 
    
    print("**************************************************", end = "\n\n")
    print(f"Region of interest:\n")
    print(f"Source name: {region_of_interest['name']}")
    print(f"Source position (ra, dec) in deg: {region_of_interest['position'].ra.deg, region_of_interest['position'].dec.deg}")
    print(f"Radius in deg: {radius_roi.value}", end = "\n\n**************************************************\n")
    
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
        pickling_catalog_roi(catalogs_roi, region_of_interest)
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


import numpy as np
# from astropy import units as u
from astropy.table import Table
from gammapy.estimators import FluxPoints
from gammapy.utils.scripts import make_path


from gammapy.datasets import FluxPointsDataset

def cut_flux_points_in_energy(datasets_counterparts, region_of_interest):
    e_ref_min = region_of_interest["e_ref_min"]
    
    e_ref_min_txt = name_to_txt(str(e_ref_min))
    source_txt  = name_to_txt(region_of_interest['name'])
    radius_roi_txt = name_to_txt(str(region_of_interest["radius_roi"]))+'degree'

    # Creates the directories to save the flux points tables 
    path_tables = get_path_tables(region_of_interest)
    
    datasets_cut_fp = Datasets()
    for index, dataset_fp in enumerate(datasets_counterparts):
        ds_name = dataset_fp.name
        flux_points= dataset_fp.data
        try:
            mask_energy = np.zeros(len(flux_points.to_table()), dtype=bool)

            for m, e_ref in enumerate(flux_points.energy_ref):
                if e_ref >= e_ref_min:
                    mask_energy[m] = True

            flux_points_mask = flux_points.to_table()[mask_energy]
            flux_points_energy = FluxPoints.from_table(flux_points_mask)

            ds = FluxPointsDataset(
                models = dataset_fp.models[0],
                data=flux_points_energy, 
                name = ds_name
            )

            datasets_cut_fp.append(ds)
            table = ds.data.to_table(
                sed_type = cfg.sed_type_e2dnde,
                formatted = True
            )    
            file_name = name_to_txt(ds_name)            
            write_tables_csv(table, path_tables, file_name)

        except Exception as error:
            print(f"The error ({dataset_fp.name}) is: {error}") 
    return datasets_cut_fp


# In[ ]:


import pickle

def pickling_catalog_roi(catalogs_roi, region_of_interest):        
    
    # Creates the directory to save the list of catalogs roi  
#     path_analysis, path_file = mkdir_sub_directory(cfg.dir_analysis, cfg.dir_catalogs_roi)
    path_file = get_path_catalogs_roi()
    
    file_name = f"catalog_{create_roi_name(region_of_interest['name'], region_of_interest['radius_roi'], region_of_interest['e_ref_min'])}"
    
    path_os = os.path.abspath(
        os.path.join(
            f"{path_file}/{file_name}{cfg.format_dat}"
        )
    )

    if path_os not in sys.path:
        sys.path.append(path_os)       

    with open(path_os, "wb") as fp:  
        pickle.dump(catalogs_roi, fp)
        
    return


# In[ ]:


import pickle
    
def unpickling_catalog_roi(region_of_interest):        
    # Creates the directory to save the catalogs roi list  
#     path_analysis, path_file = mkdir_sub_directory(cfg.dir_analysis, cfg.dir_catalogs_roi)
    path_file = get_path_catalogs_roi()
    file_name = f"catalog_{create_roi_name(region_of_interest['name'], region_of_interest['radius_roi'], region_of_interest['e_ref_min'])}"
    path_os = os.path.abspath(
        os.path.join(
            f"{path_file}/{file_name}{cfg.format_dat}"
        )
    )

    if path_os not in sys.path:
        sys.path.append(path_os)       

    with open(path_os, "rb") as fp:  
        catalogs_roi = pickle.load(fp)
    return catalogs_roi


# In[7]:





# In[12]:





# In[ ]:


def get_path_analysis():
    return mkdir_sub_directory(cfg.dir_analysis)

def get_path_catalogs_roi():
    path_analysis, path_catalogs_roi = mkdir_sub_directory(str(get_path_analysis()), cfg.dir_catalogs_roi)
    return path_catalogs_roi

def get_path_datasets(region_of_interest= None):
    if region_of_interest:
        path_analysis, path_datasets = mkdir_sub_directory(str(get_path_analysis()), f"{cfg.dir_datasets}/{region_of_interest['roi_name']}" )
    else: 
        path_analysis, path_datasets = mkdir_sub_directory(str(get_path_analysis()), f"{cfg.dir_datasets}" )

    return path_datasets

def get_path_models(region_of_interest = None):
    if region_of_interest:
        path_analysis, path_models = mkdir_sub_directory(str(get_path_analysis()), f"{cfg.dir_models}/{region_of_interest['roi_name']}")
    else:
        path_analysis, path_models = mkdir_sub_directory(str(get_path_analysis()), f"{cfg.dir_models}")        
    return path_models

def get_path_tables(region_of_interest = None):
    if region_of_interest:
        path_analysis, path_tables = mkdir_sub_directory(str(get_path_analysis()), f"{cfg.dir_tables}/{region_of_interest['roi_name']}")
    else:
        path_analysis, path_tables = mkdir_sub_directory(str(get_path_analysis()), f"{cfg.dir_tables}")
    return path_tables

def get_path_figures():
    path_analysis, path_figures = mkdir_sub_directory(str(get_path_analysis()), cfg.dir_figures)
    return path_figures

def get_path_SED_from_catalogs(region_of_interest):
    path_figures, path_SED_from_catalogs = mkdir_sub_directory(str(get_path_figures()), f"{cfg.dir_SED_from_catalogs}/{region_of_interest['roi_name']}")
    return path_SED_from_catalogs

def get_path_SED(region_of_interest):
    path_figures, path_SED = mkdir_sub_directory(str(get_path_figures()), f"{cfg.dir_SED}/{region_of_interest['roi_name']}")
    return path_SED

def get_path_flux_points(region_of_interest):
    path_figures, path_flux_points = mkdir_sub_directory(str(get_path_figures()), f"{cfg.dir_flux_points}/{region_of_interest['roi_name']}")
    return path_flux_points


# In[ ]:


from pathlib import Path


# In[ ]:


def mkdir_sub_directory(parent_directory = None, child_directory = None):
    '''Creates a directory: parent_directory/child_directory and returs the path 
    >>>mkdir_sub_directory(parent_directory, directory)
    path_parent, path_child
    '''
    if child_directory is None:

        path_parent = Path(f"{parent_directory}")
        path_parent.mkdir(exist_ok=True)
#         print("Directory '% s' created" % path_parent)
        
        return path_parent
    
    else:
        
        path_parent = Path(f"{parent_directory}")
        path_parent.mkdir(exist_ok=True)

        path_child = Path(f"{path_parent}/{child_directory}")
        path_child.mkdir(parents=True, exist_ok=True)
#         print("Directory '% s' created" % path_child)

        return (path_parent, path_child)


# In[ ]:


def name_to_txt(name):
    '''Given a `string`, `find` and `replace` the space by "_" and . by "dot"
    >>> name_to_txt(n ame.)
    namedot
    '''
    return name.replace(" ", "_").replace(".", "dot").replace(":", "")


# In[ ]:





# In[ ]:


def number_to_txt(num):
    '''Given a `number`, return a string
    >>> number_to_txt(num = 1.002222):
    1
    '''
    return ("%.2f" % num).rstrip('0').rstrip('.').replace(".", "dot")
    


# In[ ]:





# In[ ]:





# In[ ]:


from gammapy.estimators import FluxPoints

def ds_fp_from_table_fp(table, sky_model, source_name, sed_type = "e2dnde"):
    '''Returns the flux points dataset from the flux points table 
    
    >>> ds_fp_from_table_fp(table, sky_model, sed_type)
    ds_fp
    '''
    flux_points = FluxPoints.from_table(table = table, reference_model = sky_model, sed_type=sed_type)
    
    ds_name = f"{source_name}"  
    ds_fp = FluxPointsDataset(
        models = sky_model,
        data   = flux_points, 
        name   = ds_name
    )
    return ds_fp


# In[ ]:


import sys, os

def write_tables_csv(table, path_file, file_name):
# Writes the flux points table in the csv format
    path_os = os.path.abspath(
        os.path.join(
            f"{path_file}/{file_name}{cfg.format_csv}"
        )
    )

    if path_os not in sys.path:
        sys.path.append(path_os)

    table.write(
        f"{path_os}",
        format = 'ascii.ecsv', 
        overwrite = True
    )   
    return


# In[ ]:


import sys, os

def write_tables_fits(table, path_file, file_name):
    # Writes the flux points table in the fits format
    path_os = os.path.abspath(
        os.path.join(
            f"{path_file}/{file_name}{cfg.format_fits}"
        )
    )      

    if path_os not in sys.path:
        sys.path.append(path_os)

    table.write(
        f"{path_os}",
        format = 'fits', 
        overwrite = True
    )   
    return


# In[ ]:





# In[ ]:


# module_path = os.path.abspath(f'{path_my_modules}/lhaaso')
# if module_path not in sys.path:
#     print(sys.path.append(module_path))
#     sys.path.append(module_path)

# import lhaaso

