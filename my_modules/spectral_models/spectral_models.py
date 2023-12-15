from gammapy.modeling.models import (
    SkyModel, 
    PowerLawSpectralModel, 
    ExpCutoffPowerLawSpectralModel,
    LogParabolaSpectralModel,
    BrokenPowerLawSpectralModel,
)

from my_modules.utilities.utilities import name_to_txt


spectral_models_curve = {
    "PowerLawSpectralModel": PowerLawSpectralModel,
    "ExpCutoffPowerLawSpectralModel": ExpCutoffPowerLawSpectralModel,
    "LogParabolaSpectralModel": LogParabolaSpectralModel,
    "BrokenPowerLawSpectralModel": BrokenPowerLawSpectralModel,
}

class SkyModelFactory:
    def create(
            self,
            spectral_model_name, 
            datasets_names=None, 
            **kwargs,
        ):
        try: 
            spectral_model = spectral_models_curve[spectral_model_name]
            spectral_model = spectral_model(**kwargs)
        except KeyError as err:
            print(f"Model {spectral_model_name} not found")
            raise err
        except Exception as err:
            raise err
        
        if datasets_names:
            return SkyModel(
                spectral_model=spectral_model,
                datasets_names=datasets_names,
                name = f"{name_to_txt(datasets_names)}_{spectral_model.tag[1]}",
            )
        return SkyModel(spectral_model=spectral_model)
