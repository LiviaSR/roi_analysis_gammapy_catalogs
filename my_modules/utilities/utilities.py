import os
import sys
import csv
import pickle

from typing import List
from pathlib import Path

import pandas as pd 
import numpy as np


from astropy import units as u
from astropy.coordinates import SkyCoord

from gammapy.estimators import FluxPoints
from gammapy.catalog import CATALOG_REGISTRY 
from gammapy.modeling import Fit
from gammapy.datasets import (
    Datasets, 
    FluxPointsDataset,
)
from gammapy.modeling.models import (
    Models, 
    SkyModel,
)

import my_modules.config.cfg as cfg
import my_modules.cta_simulation.cta_simulation as cta

from models.system import SpiderSystem


def fit_datasets(datasets, sky_model):
    sky_model = sky_model.copy(name=sky_model.name)
    datasets.models = sky_model

    fit_result = Fit().run(datasets=datasets)

    return sky_model, fit_result


# Read a .csv file and creates a list with all systems
def read_systems_file(filename: str) -> List[SpiderSystem]:
    system_list = []
    with open(filename, 'r') as file:
        csv_reader = csv.reader(file, delimiter=";")
        next(csv_reader)  # Skip header line
        for row in csv_reader:
            system_list.append(SpiderSystem(
                source_name=row[0],
                pos_ra=float(row[4]),
                pos_dec=float(row[5]),
            ))

    return system_list

# Create empty directories of each systems inside the
#  <base_dir> directory. 
def create_output_dir(base_dir, name):
    output_dir = f"{base_dir}/{name}"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        
    return output_dir


def set_source_info(system_data):
    source_name = system_data.source_name
    pos_ra = system_data.pos_ra*u.deg
    pos_dec = system_data.pos_dec*u.deg

#     if any([source_RA.unit !=  cfg.unit_deg, source_dec.unit !=  cfg.unit_deg]):
#         raise Exception("Sorry, there is a error: celestial coordinates (RA, dec.) units is not in degrees") 
    
    return  {
        'name': source_name,
        'position': SkyCoord(pos_ra, pos_dec) 
    }



def write_datasets_models(datasets,region_of_interest, directory_name, path_datasets = None):   
    if path_datasets == None:
        path_datasets = get_path_datasets(region_of_interest)        
    
    path_datasets, path_file = mkdir_sub_directory(str(path_datasets), directory_name)
    datasets.write(filename=f"{path_file}/datasets{cfg.format_yaml}", filename_models=f"{path_file}/models{cfg.format_yaml}", overwrite=True)
    return


def read_datasets_models(region_of_interest, directory_name, path_datasets = None):
    if path_datasets == None:
        path_datasets = get_path_datasets(region_of_interest)  
    path_datasets, path_file = mkdir_sub_directory(str(path_datasets), directory_name)
    return Datasets.read(filename=f"{path_file}/datasets{cfg.format_yaml}", filename_models=f"{path_file}/models{cfg.format_yaml}")



def create_file_name(pulsar_info, irf_name, skymodel_name): 
    """
    """
    return f"{name_to_txt(pulsar_info['name'])}_irf_{irf_name}_{skymodel_name}"


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

def create_region_of_interest(
        source_info=None, 
        radius_roi=1, 
        e_ref_min=None, 
        e_ref_max=None
    ):
    radius_roi = radius_roi * u.deg

    region_of_interest = {
        **source_info,
        "radius_roi": radius_roi,
        "e_ref_min": e_ref_min,
        "e_ref_max": e_ref_max,
        "roi_name": create_roi_name(source_info['name'], radius_roi, e_ref_min, e_ref_max),
    }
    
    return region_of_interest


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


def get_dict_sources_gammapy(sources_gammapy, datasets_gammapy):
    dict_sources_gammapy = {}
    for index, (source, dataset) in enumerate(zip(sources_gammapy, datasets_gammapy)):
        dict_sources_gammapy[dataset.name] = {'position': source.position}
    return dict_sources_gammapy


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
            if counterpart.flux_points_table is None or counterpart.spectral_model() is None:
                continue
            n_counterparts+=1   
            counterpart_name = counterpart.name            
            counterpart_flux_points = counterpart.flux_points
            n_flux_points+=1
            counterpart_spectral_model = counterpart.spectral_model()
    
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
                
                            
    print(f"Total number of counterparts: {n_counterparts}")
    print(f"Total number of flux points tables: {n_flux_points}")
    return counterparts, datasets_counterparts, models_counterparts


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

    with open(path_os, "rb") as fp:  
        catalogs_roi = pickle.load(fp)
    return catalogs_roi


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


def name_to_txt(name):
    '''Given a `string`, `find` and `replace` the space by "_" and . by "dot"
    >>> name_to_txt(n ame.)
    namedot
    '''
    return name.replace(" ", "_").replace(".", "dot").replace(":", "")


def number_to_txt(num):
    '''Given a `number`, return a string
    >>> number_to_txt(num = 1.002222):
    1
    '''
    return ("%.2f" % num).rstrip('0').rstrip('.').replace(".", "dot")

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
