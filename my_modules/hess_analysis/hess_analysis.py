from gammapy.catalog import CATALOG_REGISTRY 
from gammapy.modeling.models import SkyModel

import my_modules.config.cfg as cfg
import my_modules.utilities.utilities as utl
import my_modules.spectral_models.spectral_models as spec

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

