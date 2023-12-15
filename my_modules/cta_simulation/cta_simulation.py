# astropy affiliated packages imports
from regions import CircleSkyRegion
from astropy import units as u

from gammapy.irf import load_cta_irfs
from gammapy.maps import MapAxis, RegionGeom
from gammapy.modeling import Fit
from gammapy.estimators import FluxPointsEstimator
from gammapy.stats import WStatCountsStatistic
from gammapy.datasets import (
    Datasets, 
    SpectrumDataset,
    SpectrumDatasetOnOff,
)
from astropy.coordinates import (
    SkyCoord, 
    Distance,
    Angle,
)
from gammapy.data import (
    Observation,
    observatory_locations,
)
from gammapy.makers import (
    SafeMaskMaker,
    SpectrumDatasetMaker, 
)
from gammapy.modeling.models import (
    Models, 
    SkyModel,
    EBLAbsorptionNormSpectralModel,
)

import my_modules.config.cfg as cfg

def get_dict_pulsars():  
    dict_pulsars = {
        'PSR J1826-1334': {
            'position': SkyCoord(276.554896, -13.57967, unit="deg"), 
            'distance': 3.1 * u.kpc,
            'age': 21.4 * u.kyr,
            'luminosity': 2.8e+36 * u.Unit("erg s-1")
        },
         'PSR J1826-1256': {
            'position': SkyCoord(276.53554, -12.94250, unit="deg"), 
            'distance': 1.6 * u.kpc,
             'age': 14.4 * u.kyr,
             'luminosity': 3.6e+36 * u.Unit("erg s-1")
        },
         'PSR J1837-0604': {
            'position': SkyCoord(279.43146, -6.0803, unit="deg"), 
            'distance': 4.8 * u.kpc,
            'age': 33.8 * u.kyr,
            'luminosity': 2.0e+36 * u.Unit("erg s-1")
        },
         'PSR J1838-0537': {
            'position': SkyCoord(279.73342, -5.6192, unit="deg"), 
            'distance': 1.3 * u.kpc,
             'age': 4.9 * u.kyr,
             'luminosity': 6.0e+36 * u.Unit("erg s-1")
        },
    }
    return dict_pulsars


def set_pulsar_info(dict_pulsars, pulsar_index):
    """
    Sets the pulsar info into a dictionary
    """
    pulsar_name = list(dict_pulsars.keys())[pulsar_index]
    position_RA = list(dict_pulsars.values())[pulsar_index]["position"].ra
    position_dec = list(dict_pulsars.values())[pulsar_index]["position"].dec
    pulsar_pos = SkyCoord(position_RA, position_dec) # Source Position
    pulsar_dist = Distance(list(dict_pulsars.values())[pulsar_index]["distance"])
    pulsar_red = float(pulsar_dist.compute_z()) # The source redshift for this distance assuming its physical distance is a luminosity distance.

    return  {
        "name": pulsar_name,
        "position": pulsar_pos,
        "distance": pulsar_dist,
        "redshift": pulsar_red
    }

def get_pointing(pulsar_info, livetime, offset):
    return SkyCoord(pulsar_info['position'].ra, pulsar_info['position'].dec + offset, unit="deg")


# In this simulation, we use the CTA-1DC irfs shipped with gammapy
path_irfs_files = '/home/born-again/Documents/GitHub/gammapy/gammapy-notebooks/0.20.1/tutorials/data/caldb/data/cta/prod3b-v2/bcf'

irf_z = [20,40,60]
irf_h = [0.5, 5, 50]
irf_loc = [("cta_north", "North"),("cta_south", "South")]


def create_irf_name(irf_zenith = 0, irf_hours = 1, irf_site = 1):
    return f'{irf_loc[irf_site][1]}_z{irf_z[irf_zenith]}_{irf_h[irf_hours]}h'


def load_irfs(irf_name):
    irf_filename = f"{path_irfs_files}/{irf_name}/irf_file{cfg.format_fits}"
    return load_cta_irfs(irf_filename)


def create_location(irf_site):
    return observatory_locations[irf_loc[irf_site][0]]


def create_observation(location,pointing,livetime,irfs):
    observation = Observation.create(
        pointing=pointing,
        livetime=livetime,
        irfs=irfs,
        location=location,
    )
    return observation


def create_path_model(region_of_interest):
    return f"/home/born-again/Documents/GitHub/CTA_projects/my_notebooks/counterparts/analysis/models/{region_of_interest['roi_name']}"



def get_spectral_model(region_of_interest, model_name):
    path_model = create_path_model(region_of_interest)
    models = Models.read(f"{path_model}/{model_name}{cfg.format_yaml}")
    return models[0].spectral_model


def get_absorption_model(pulsar_info, abs_model = 'franceschini'):
    # Available models in gammapy-data:{'franceschini', 'dominguez', 'finke'}
    absorption = EBLAbsorptionNormSpectralModel.read_builtin(
        reference = abs_model, 
        redshift=pulsar_info['redshift']
    )
    return absorption

    
def create_abs_sky_model(spectral_model, absorption, model_name):
    absspecmodel =  spectral_model * absorption # CompoundSpectralModel
    sky_model=SkyModel(
#     spectral_model=spectral_model, 
        spectral_model=absspecmodel, 
        name=model_name
    )
    return sky_model

 
def set_reconstructed_energy_axis(e_edges_min, e_edges_max):    
    # Reconstructed energy axis
    energy_reco = MapAxis.from_energy_bounds(
        e_edges_min, 
        e_edges_max, 
        nbin=5, 
        per_decade=True, 
        name="energy"
    )
    return energy_reco


def define_on_region(pulsar_info, radius):
    on_region_radius = Angle(radius)
    on_region = CircleSkyRegion(
        center=pulsar_info["position"], 
        radius=on_region_radius
    )
    return on_region


def define_geometry(on_region, energy_reco):
    #Defines the geometry:
    geom = RegionGeom.create(
        region=on_region, 
        axes=[energy_reco]
    )
    return geom


def set_true_energy_axis(e_edges_min, e_edges_max):
    # Defines the true energy axis:
    # true energy axis should be wider than reco energy axis
    energy_true = MapAxis.from_energy_bounds(
        0.3*e_edges_min, 
        3*e_edges_max, 
        nbin=8, 
        per_decade=True, 
        name="energy_true"
    )
    return energy_true


def spectrum_dataset_maker():    
    # Make spectrum for a single IACT observation:
    # The irfs and background are computed at a single fixed offset, which is recommended only for point-sources.
    maker = SpectrumDatasetMaker(
    #     containment_correction=True, # Apply containment correction for point sources and circular on regions.
        selection=["edisp", "background", "exposure"], # Selecting which maps to make
        use_region_center=False
    )
    return maker

def safe_mask_maker():
    # Make safe data range mask for a given observation.
    return SafeMaskMaker(methods=["bkg-peak"], aeff_percent=10) 


def create_dataset_empty(geom, energy_true):
    # Create a MapDataset object with zero filled maps.
    dataset_empty = SpectrumDataset.create(
        geom=geom, 
        energy_axis_true=energy_true,
        name="obs-0"
    )
    return dataset_empty


def make_map_dataset(maker, safe_maker,dataset_empty, observation,skymodel):
    # Make map dataset:
    dataset = maker.run(dataset_empty, observation) 
    dataset = safe_maker.run(dataset, observation)

    # Set the model on the dataset, and fake
    dataset.models = skymodel
    dataset.fake(random_state=42)
    return dataset


def spectrum_dataset_on_off(dataset):
    # Spectrum dataset for on-off likelihood fitting.
    dataset_onoff = SpectrumDatasetOnOff.from_spectrum_dataset(
        dataset=dataset, 
        acceptance=1, 
        acceptance_off=5
    )

    # Simulate fake counts (on and off) for the current model and reduced IRFs.
    dataset_onoff.fake(
        random_state='random-seed', 
        npred_background=dataset.npred_background()
    )
    return dataset_onoff


def compute_significance(dataset_onoff):
    # Class to compute statistics for Poisson distributed variable with unknown background.
    significance = WStatCountsStatistic(
        n_on=sum(dataset_onoff.counts.data), 
        n_off=sum(dataset_onoff.counts_off.data), 
        alpha=0.2).sqrt_ts
    return significance


def simulate_observations(n_obs,dataset, dataset_onoff):    
    datasets = Datasets()

    for idx in range(n_obs):
        dataset_onoff.fake(
            random_state=idx, 
            npred_background=dataset.npred_background()
        )
        dataset_fake = dataset_onoff.copy(name=f"obs-{idx}")
        dataset_fake.meta_table["OBS_ID"] = [idx]
        datasets.append(dataset_fake)
    return datasets


def compute_flux_points(datasets, skymodel, e_edges_min,e_edges_max):
    #Compute flux points
    datasets.models = [skymodel]

    fit_joint = Fit()
    result_joint = fit_joint.run(datasets=datasets)

    # we make a copy here to compare it later
    model_best_joint = skymodel.copy()

    energy_edges = MapAxis.from_energy_bounds(e_edges_min,e_edges_max, nbin=12).edges

    fpe = FluxPointsEstimator(energy_edges=energy_edges, source=skymodel.name, selection_optional="all")
    flux_points = fpe.run(datasets=datasets)
    table = flux_points.to_table(sed_type="e2dnde", formatted=True)
    CTA_table = set_CTA_table_to_erg(table)

    return datasets, model_best_joint, CTA_table


def set_CTA_table_to_erg(table):
    table["e2dnde"][:] = table["e2dnde"].to(u.Unit("erg cm-2 s-1")) 
    table["e2dnde"].unit = u.Unit("erg cm-2 s-1")

    table["e2dnde_err"][:] = table["e2dnde_err"].to(u.Unit("erg cm-2 s-1")) 
    table["e2dnde_err"].unit = u.Unit("erg cm-2 s-1")

    table["e2dnde_errp"][:] = table["e2dnde_errp"].to(u.Unit("erg cm-2 s-1")) 
    table["e2dnde_errp"].unit = u.Unit("erg cm-2 s-1")

    table["e2dnde_errn"][:] = table["e2dnde_errn"].to(u.Unit("erg cm-2 s-1")) 
    table["e2dnde_errn"].unit = u.Unit("erg cm-2 s-1")

    table["e2dnde_ul"][:] = table["e2dnde_ul"].to(u.Unit("erg cm-2 s-1")) 
    table["e2dnde_ul"].unit = u.Unit("erg cm-2 s-1")
    return table

