#import pickle

import numpy as np
from pprint import pprint

from gammapy.estimators import FluxPoints
from gammapy.modeling.models import SkyModel, Models
from gammapy.datasets import (
    FluxPointsDataset,
    Datasets,
)
from gammapy.catalog.hawc import SourceCatalogObjectHWCBase


import my_modules.utilities.utilities as utl
import my_modules.hess_analysis.hess_analysis as hess

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


def cut_flux_points_energy_range(flux_points, e_ref_min, e_ref_max):
    # Raise exception of no energy cut is defined 
    if all([e_ref_min ==  None, e_ref_max ==  None]):
            raise Exception(f"Sorry, there is a error: e_ref_min is None and e_ref_max is None)") 
    
    # Create an array of zeros with a desired boolean type
    mask_energy = np.zeros(len(flux_points.to_table()), dtype=bool)

    for m, e_ref in enumerate(flux_points.energy_ref):
        # If the energy of reference is bigger than the lower energy cut, 
        # set a True value in the array
        if e_ref_min != None:
            if e_ref >= e_ref_min:
                mask_energy[m] = True
        # If the energy of reference is smaller than the upper energy cut, 
        # set a True value in the array   
        if e_ref_max != None:
            if e_ref <= e_ref_max:
                mask_energy[m] = True
    # Define a flux points table with energy values that satisfies the energy cuts
    flux_points_mask = flux_points.to_table()[mask_energy]
    if len(flux_points_mask) == 0:
        return None

    # Read the flux points from the table created above
    flux_points = FluxPoints.from_table(flux_points_mask)     
        
    return flux_points


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

        # for obj in catalog[mask_roi]:
        #     print('printing obj position')
        #     print(obj.position)
        #     print(source_position.ra.degree, source_position.dec.degree)
        #     print(source_position.separation(obj.position).degree)
        
        if len(catalog[mask_roi].table):
            catalogs_roi.append(catalog[mask_roi])
            numbers_catalogs_roi += 1
        else:
            catalogs_no_counterparts.append(f"{catalog.tag}: {catalog.description}")
    
            
    if numbers_catalogs_roi:
        utl.pickling_catalog_roi(catalogs_roi, region_of_interest)
        #print(f"\n{numbers_catalogs_roi} catalogs with sources within the region of interest:", end = "\n\n")
        for catalog in catalogs_roi:
            continue
            #print(f"{catalog.tag}: {catalog.description}")
            #display(catalog.table)
    else:
        print("No catalogs with sources in the region of interest!", end = "\n\n")

    if numbers_catalogs_roi and len(catalogs_no_counterparts):
        print("Catalogs without sources within the region of interest:", end = "\n\n")
        for index, catalog_no_counterpart in enumerate(catalogs_no_counterparts):                            
            print(catalog_no_counterpart)

    return catalogs_roi


def get_datasets_flux_points_gammapy(region_of_interest):
    '''
    Select a catalog subset (only sources within a region of interest)
    '''  
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

    for catalog in catalogs_roi:
        cat_tag = catalog.tag
        for counterpart in catalog:
            if isinstance(counterpart, SourceCatalogObjectHWCBase):
                #print("Found HWC Catalog entry:")
                #print(type(counterpart))
                #pprint(counterpart.data)
                #hawc_counterpart_flux_points = FluxPoints()
                continue
            if counterpart.flux_points_table is None or counterpart.spectral_model() is None:
                continue
            n_counterparts+=1
            flux_points_table_df = counterpart.flux_points_table.to_pandas()
            #print(counterpart.name)
            #print(flux_points_table_df)  
            #pprint(counterpart.spectral_model())
            #print(counterpart.flux_points)
            counterpart_name = counterpart.name            
            flux_points = counterpart.flux_points

            if e_ref_min != None or e_ref_max != None:
                flux_points = cut_flux_points_energy_range(flux_points, e_ref_min, e_ref_max)
                if flux_points == None:
                    print(counterpart_name)
                    print("Warning: No counterparts in the selected energy range.")
                    continue

            counterpart_spectral_model = counterpart.spectral_model()
            
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
            
            n_flux_points+=1
            models_counterparts.append(counterpart_model)  # Add the counterpart_model to models()
            
            counterparts.append(counterpart)
            datasets_counterparts.append(ds)
    
    datasets_counterparts.models = models_counterparts
            
    print('---------------------------')        
    pprint(f"Total number of counterparts: {n_counterparts}")
    pprint(f"Total number of flux points tables: {n_flux_points}")
    return counterparts, datasets_counterparts, models_counterparts

