from pathlib import Path

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord

from gammapy.modeling.models import SkyModel, Models
from gammapy.datasets import Datasets

import my_modules.config.cfg as cfg
import my_modules.utilities.utilities as utl
import my_modules.spectral_models.spectral_models as spec

def get_dict():
    '''
    Dictionary of the LHAASO PeVatrons (dict keys: Source name based on J2000 coordinates) information: 
"position": Right ascension (in degrees), Declination (in degrees) 
"flux": Differential photon fluxes at 100 TeV (in Crab Units (CU)), The corresponding error
obs: CU is the flux of the Crab Nebula at 100 TeV; 1 CU = 6.1 x 10-17 photons TeV-1 cm-2 s-1
see https://www.nature.com/articles/s41586-021-03498-z
    '''
    CU = 6.1e-17 * u.Unit("TeV-1 cm-2 s-1")
    return {
    "LHAASO J0534+2202": {
        "position": SkyCoord(83.55, 22.05, unit="deg"),
        "flux":(1.00, 0.14)* CU
    },
    "LHAASO J1825-1326": {
        "position": SkyCoord(276.45, -13.45, unit="deg"),
        "flux": (3.57,0.52)* CU
    },    
    "LHAASO J1839-0545": {
        "position": SkyCoord(279.95, -5.75, unit="deg"),
        "flux": (0.70,0.18)* CU
    },    
    "LHAASO J1843-0338": {
        "position": SkyCoord(280.75, -3.65, unit="deg"),
        "flux": (0.73,0.17)* CU
    }, 
    "LHAASO J1849-0003": {
        "position": SkyCoord(282.35, -0.05, unit="deg"),
        "flux": (0.74,0.15)* CU
    }, 
    "LHAASO J1908+0621": {
        "position": SkyCoord(287.05, 6.35, unit="deg"),
        "flux": (1.36,0.18)* CU
    },
    "LHAASO J1929+1745": {
        "position": SkyCoord(292.25, 17.75, unit="deg"),
        "flux": (0.38,0.09)* CU
    },    
    "LHAASO J1956+2845": {
        "position": SkyCoord(299.05, 28.75, unit="deg"),
        "flux": (0.41,0.09)* CU
    },
    "LHAASO J2018+3651": {
        "position": SkyCoord(304.75, 36.85, unit="deg"),
        "flux": (0.50,0.10)* CU
    },
    "LHAASO J2032+4102": {
        "position": SkyCoord(308.05, 41.05, unit="deg"),
        "flux": (0.54,0.10)* CU
    },
    "LHAASO J2108+5157": {
        "position": SkyCoord(317.15, 51.95, unit="deg"),
        "flux": (0.38,0.09)* CU
    },
    "LHAASO J2226+6057": {
        "position": SkyCoord(336.75, 60.95, unit="deg"),
        "flux": (1.05,0.16)* CU
    }
}
    
def get_dataset(source_name): 
    
    # sky model LHAASO J1825-1326
    if source_name == 'LHAASO J1825-1326':
        sky_model = spec.sky_model_lp(
            alpha = 0.92,
            amplitude = "1e-12 cm-2 s-1 TeV-1",
            reference = 10 * u.TeV,
            beta = 1.19,
            datasets_names = source_name
        )
        
    # sky model LHAASO J1908+0621
    elif source_name == "LHAASO J1908+0621":
        sky_model = spec.sky_model_lp(
        alpha = 2.27,
        amplitude = "1e-12 cm-2 s-1 TeV-1",
        reference = 10 * u.TeV,
        beta = 0.46,
        datasets_names = source_name
    )

    # sky model LHAASO J2226+6057
    elif source_name == "LHAASO J2226+6057":
        sky_model = spec.sky_model_lp(
        alpha = 1.56,
        amplitude = "1e-12 cm-2 s-1 TeV-1",
        reference = 10 * u.TeV,
        beta = 0.88,
        datasets_names = source_name
    )
    else:
        sky_model = spec.sky_model_pl(
            datasets_names = source_name
        )
     
    table = table_to_SED_format(cfg.path_fp_LHAASO, utl.name_to_txt(source_name))
        
    return utl.ds_fp_from_table_fp(table = table, sky_model = sky_model, source_name = source_name)


def table_to_SED_format(path, file_name):
    '''
    Normalization Representation
    The SED format is a flexible specification for representing one-dimensional spectra 
    (distributions of amplitude vs. energy).
    
    '''
    
    format_dat = '.dat'
    file_path = Path(f'{path}/{file_name}{format_dat}') 

    table = Table.read(file_path,format='ascii', delimiter=' ', comment='#')
    
#     display(table)

    table['col1'] = table['col1']/1e12
    table.rename_column('col1', 'e_ref')
    table['e_ref'].unit = u.TeV

    #     table['col5'] = table['col5']/1e12
    #     table.rename_column('col5', 'e_min')
    #     table['e_min'].unit = u.TeV

    #     table['col6'] = table['col6']/1e12
    #     table.rename_column('col6', 'e_max')
    #     table['e_max'].unit = u.TeV
    
    table['col2'] = table['col2']*u.erg.to("TeV")


    table.rename_column('col2', 'e2dnde')
    table['e2dnde'].unit = u.Unit("TeV cm-2 s-1")
    
    table['col3'] = table['col3']*u.erg.to("TeV")
    table.rename_column('col3', 'e2dnde_errp')
    table['e2dnde_errp'].unit = u.Unit("TeV cm-2 s-1")
    
    table['col4'] = table['col4']*u.erg.to("TeV")
    table.rename_column('col4', 'e2dnde_errn')
    table['e2dnde_errn'].unit = u.Unit("TeV cm-2 s-1")

    table.meta["SED_TYPE"] = "e2dnde"
    table.meta["name"] = "table"
    try:
        table.remove_columns(['col5', 'col6'])
    except:
        pass
    
    display(table)
    
    return table


def get_LHAASO_tables_datasets(dict_LHAASO):  
    # sky model LHAASO J1825-1326
    source_name = 'LHAASO J1825-1326'
    sky_model_LHAASO_J1825 = spec.sky_model_lp(
        alpha = 0.92,
        amplitude = "1e-12 cm-2 s-1 TeV-1",
        reference = 10 * u.TeV,
        beta = 1.19,
        datasets_names = source_name
    )
    # sky_model_LHAASO_J1825.spectral_model.plot(energy_bounds = [3e1, 3e3] * u.TeV, sed_type=cfg.sed_type_e2dnde)

    # sky model LHAASO J1908+0621
    source_name = "LHAASO J1908+0621"
    sky_model_LHAASO_J1908 = spec.sky_model_lp(
        alpha = 2.27,
        amplitude = "1e-12 cm-2 s-1 TeV-1",
        reference = 10 * u.TeV,
        beta = 0.46,
        datasets_names = source_name
    )

    # sky model LHAASO J2226+6057
    source_name = "LHAASO J2226+6057"
    sky_model_LHAASO_J2226 = spec.sky_model_lp(
        alpha = 1.56,
        amplitude = "1e-12 cm-2 s-1 TeV-1",
        reference = 10 * u.TeV,
        beta = 0.88,
        datasets_names = source_name
    )
    
    tables = []
    datasets = []
    models = Models()
    datasets = Datasets()
    for source_index, source_name in enumerate(list(dict_LHAASO.keys())):
        print(source_index, source_name)
        table = table_to_SED_format(cfg.path_fp_LHAASO, utl.name_to_txt(source_name))
        tables.append(table)
        if source_name == 'LHAASO J1825-1326':
            sky_model = sky_model_LHAASO_J1825.copy(name = sky_model_LHAASO_J1825.name, datasets_names = source_name)
        elif source_name == 'LHAASO J1908+0621':
            sky_model = sky_model_LHAASO_J1908.copy(name = sky_model_LHAASO_J1908.name, datasets_names = source_name)
        elif source_name == 'LHAASO J2226+6057':
            sky_model = sky_model_LHAASO_J2226.copy(name = sky_model_LHAASO_J2226.name, datasets_names = source_name)
        else:
            sky_model = SkyModel(
                spectral_model = spec.sky_model_pl().spectral_model, 
                name = f"{utl.name_to_txt(source_name)}_pl",
                datasets_names = source_name
            )    
        models.append(sky_model)
        dataset = utl.ds_fp_from_table_fp(table = table, sky_model = sky_model, source_name = source_name)
        datasets.append(dataset)
        
    return tables, datasets, models
