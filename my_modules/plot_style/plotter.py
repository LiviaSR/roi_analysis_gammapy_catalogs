#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from gammapy.modeling.models import Models
from gammapy.datasets import Datasets


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


# In[4]:


module_path = os.path.abspath(f'{path_my_modules}/{cfg.dir_utilities}')
if module_path not in sys.path:
    sys.path.append(module_path)

import utilities as utl
importlib.reload(utl)


# In[ ]:


def savefig(path_file, file_name):
    ''' Saves figures (.pdf, .png  and .svg) in the path_child directoty    
    savefig(path_child, child_name)
    >>> plt.savefig(file, bbox_inches='tight')
    '''
    formats_file = [cfg.format_pdf, cfg.format_png, cfg.format_svg]
    for format_file in formats_file: 
        file = path_file / f'{file_name}{format_file}'
        plt.savefig(file, bbox_inches='tight')


# In[ ]:


def set_leg_style(dict_leg_style, datasets = None, models = None, color = None, marker = None, linestyle = None):
    if all([datasets ==  None, models ==  None]):
        return print("Sorry, there is error: 'datasets =  None' and 'models =  None'")
    else: 
#             marker_ds = marker
#             color_ds = color
    
        if datasets !=  None:
            dict_leg_style = set_leg_style_datasets(dict_leg_style, datasets, color, marker)
        
        if models !=  None:
            dict_leg_style = set_leg_style_models(dict_leg_style, models, color, linestyle)
        
        return dict_leg_style


# In[ ]:


def set_leg_style_datasets(dict_leg_style, datasets, color = None, marker = None):
    datasets = Datasets(datasets)
    marker_ds = marker
    color_ds = color
    if not marker_ds:
        while len(cfg.markers) < len(datasets) +1:
            cfg.markers.extend(cfg.markers)
    if not color_ds:      
        while len(cfg.colors) < len(datasets) +1:
            cfg.colors.extend(cfg.colors)

    for index, dataset in enumerate(datasets):
        if not color_ds:
            color = cfg.colors[index]

        if not color_ds:
            marker = cfg.markers[index]
        
        #############################
        if dataset.name.find('LHAASO') != -1:
            color = cfg.color_lhaaso
            marker = cfg.marker_lhaaso
            
        if dataset.name.find('CTA') != -1:
            color = cfg.color_cta
            marker = cfg.marker_cta
        #############################    
        dict_leg_style[dataset.name] = (color, marker)
    return dict_leg_style


# In[ ]:


def set_leg_style_models(dict_leg_style, models, color = None, linestyle = None):
    models = Models(models)
    color_m = color
    linestyle_m = linestyle
    
    if not linestyle:
        while len(cfg.linestyles) < len(models) +1:
            cfg.linestyles.extend(cfg.linestyles)
    if not color_m:      
        while len(cfg.colors) < len(models) +1:
            cfg.colors.extend(cfg.colors)

    for index, model in enumerate(models):
        if not color_m:
            color = "black"
            
        linestyle = cfg.linestyles[index]
        dict_leg_style[model.name] = (color, linestyle)
    return dict_leg_style


# In[ ]:


import re
def set_label_datasets(dataset_name):
    
    test_string = dataset_name
    spl_word = ':'
    match = re.search(spl_word, test_string)
    if match:
        source_name = test_string[:match.end()-1]
        cat_tag = (test_string[match.end():]).replace(" ", "")
        if cat_tag == 'hgps':
            year = "2018"
        else:
            if source_name == "HESS J1825-137":
                year = "2006"
            if source_name == "HESS J1826-130":
                year = "2017"
            if source_name == "HESS J1837-069":
                year = "2006"
            if source_name == "HESS J1841-055":
                year = "2018b"
        label = f'{source_name} ({year})' 
    else:
        label = dataset_name

    return label


# In[ ]:


from astropy import units as u
def plot_SED(
    name = "region_of_interest", 
    datasets = None,  
    models = None,
    dict_leg_style = None, 
    region_of_interest = None,
    sed_type = "e2dnde", 
    dict_plot_axis =  dict(
    label =  (r'$\rm{E\ [TeV] }$', r'$\rm{E^2\ J(E)\ [TeV\ cm^{-2}\ s^{-1}] }$'),
    units =  (          'TeV',                       'TeV  cm-2     s-1')
),
    dict_plot_limits = dict(
        energy_bounds = [1e-5, 3e2] * u.TeV,
        ylim = [1e-23, 1e-7]
    ),
    dict_leg_place = dict(
#         bbox_to_anchor = (0, -0.45), # Set legend outside plot
        ncol=3, 
        loc='lower left', 
    )
):    
    
    ax = plt.subplot()
    
    ax.xaxis.set_units(u.Unit(dict_plot_axis['units'][0]))
    ax.yaxis.set_units(u.Unit(dict_plot_axis['units'][1]))

    kwargs = {
        "ax": ax, 
        "sed_type": sed_type,
#         "uplims": True
    }

#     while len(cfg.markers) < len(datasets) + 1:
#          cfg.markers.extend( cfg.markers)
                        
    for index, dataset in enumerate(datasets):
        color = dict_leg_style[dataset.name][0]
        marker = dict_leg_style[dataset.name][1]
        
        
        
        label = set_label_datasets(dataset.name) 
            
        dataset.data.plot(
                    label = label, 
                    marker = marker, 
                    color=color,
                    **kwargs
                )
    if models:
        path_file =  utl.get_path_SED(region_of_interest)  
        for index, model in enumerate(models):
            linestyle = dict_leg_style[model.name][1]
            color = dict_leg_style[model.name][0]
            spectral_model = model.spectral_model
            
#             spectral_model.plot(label = f"{model.name} (fit)", energy_bounds=dict_plot_limits['energy_bounds'],   marker = ',', color="black", **kwargs)
            energy_bounds = [7e-2, 8e2] * u.TeV
#             energy_bounds=dict_plot_limits['energy_bounds']
            spectral_model.plot(energy_bounds=energy_bounds,  linestyle = linestyle, marker = ',', color=color, **kwargs)

            spectral_model.plot_error(energy_bounds=energy_bounds,**kwargs)
    else:
        path_file =  utl.get_path_flux_points(region_of_interest)  
    ax.set_ylim(dict_plot_limits['ylim'])
    ax.set_xlim(dict_plot_limits['energy_bounds'])
    
    ax.legend(**dict_leg_place)
    
    plt.xlabel(dict_plot_axis['label'][0])   
    plt.ylabel(dict_plot_axis['label'][1])
    
    file_name = utl.name_to_txt(name)
        
    savefig(path_file, file_name)
    
#     plt.savefig(path_file, bbox_inches='tight')
#     plt.grid(which="both")
    plt.show()
    
    return


# In[5]:


from astropy import units as u
import matplotlib.pyplot as plt # A collection of command style functions

def SED_from_catalogs(
    counterparts, datasets_counterparts, models_counterparts, colors_dict, region_of_interest,
                           sed_type = "e2dnde", 
                           axis_dict =  dict(
    label =  (r'$\rm{E\ [TeV] }$', r'$\rm{E^2\ J(E)\ [TeV\ cm^{-2}\ s^{-1}] }$'),
    units =  (          'TeV',                       'TeV  cm-2     s-1')
),                         
    energy_bounds = [1e-5, 1e2] * u.TeV, 
    ylim = [1e-13, 1e-9]
                          ):
    
    sed_type = sed_type
    axis_dict =axis_dict
    energy_bounds =energy_bounds
    ylim = ylim
    
        
    for index, (counterpart, dataset, model) in enumerate(zip(counterparts, datasets_counterparts, models_counterparts)):    

        counterpart_name = counterpart.name
        flux_points = counterpart.flux_points
#         spectral_model = model.spectral_model
        spectral_model = counterpart.spectral_model()
        spectral_model_tag = spectral_model.tag[0]
        spectral_model_tag_short = spectral_model.tag[1]

        ax = plt.subplot()
        ax.xaxis.set_units(u.Unit(axis_dict['units'][0]))
        ax.yaxis.set_units(u.Unit(axis_dict['units'][1]))
    
        xlabel = axis_dict['label'][0]
        ylabel = axis_dict['label'][1]

        kwargs = {
            "ax": ax, 
            "sed_type": sed_type
        }
        kwargs_fit = {
            "label": f"{spectral_model_tag_short} (fit)"
        }

        color = colors_dict[dataset.name][0]
        marker = colors_dict[dataset.name][1]
        flux_points.plot(label = dataset.name, color= color, marker = marker, **kwargs)

        energy_bounds = flux_points.energy_min[0], flux_points.energy_max[-1]
    #     e2dnde_errn = flux_points.e2dnde_errn.data
    #     e2dnde_errp = flux_points.e2dnde_errp.data
    #     ylim =min(np.nanmin(e2dnde_errp),np.nanmin(e2dnde_errn)), max(np.nanmax(e2dnde_errp),np.nanmax(e2dnde_errn))

        spectral_model.plot(energy_bounds=energy_bounds, marker = ',', ls= "-", color="k", **kwargs, **kwargs_fit)
        spectral_model.plot_error(energy_bounds=energy_bounds, **kwargs)

    #     kwargs_spectrum = {"kwargs_model": {"color":"red", "ls":"--"}, "kwargs_fp":{"color":"green", "marker":"o"}}
    #     dataset_fp.plot_spectrum(**kwargs, **kwargs_spectrum)  
    #    ax.set_ylim(ylim)
    
        ax.set_xlim(energy_bounds)
        plt.xlabel(xlabel)   
        plt.ylabel(ylabel)
        plt.legend()
        file_name = f'{utl.name_to_txt(dataset.name.replace(":",""))}_{spectral_model_tag_short}{cfg.format_png}'
        file_path = utl.get_path_SED_from_catalogs(region_of_interest) / file_name 
        plt.savefig(file_path, bbox_inches='tight')
        plt.show()
      


# In[ ]:


from astropy import units as u
def SED_flux_points(
    datasets = None,  
    spectral_model = None,
    colors_dict = None, 
    region_of_interest = None,
    sed_type = "e2dnde", 
    axis_dict =  dict(
    label =  (r'$\rm{E\ [TeV] }$', r'$\rm{E^2\ J(E)\ [TeV\ cm^{-2}\ s^{-1}] }$'),
    units =  (          'TeV',                       'TeV  cm-2     s-1')
),
    limits_dict = dict(
        energy_bounds = [1e-5, 3e2] * u.TeV,
        ylim = [1e-23, 1e-7]
    ),
    legend_dict = dict(
#         bbox_to_anchor = (0, -0.45), # Set legend outside plot
        ncol=3, 
        loc='lower left', 
    )
):    
    
    ax = plt.subplot()
    
    ax.xaxis.set_units(u.Unit(axis_dict['units'][0]))
    ax.yaxis.set_units(u.Unit(axis_dict['units'][1]))

    kwargs = {
        "ax": ax, 
        "sed_type": sed_type,
#         "uplims": True
    }

#     while len(cfg.markers) < len(datasets) + 1:
#          cfg.markers.extend( cfg.markers)
                        
    for index, dataset in enumerate(datasets):
        color = colors_dict[dataset.name][0]
        marker = colors_dict[dataset.name][1]
        dataset.data.plot(
                    label = dataset.name, 
                    marker = marker, 
                    color=color,
                    **kwargs
                )
    if spectral_model:
        spectral_model.plot(label = "Counterparts (fit)", energy_bounds=limits_dict['energy_bounds'],  marker = ',', color="k", **kwargs)
        spectral_model.plot_error(energy_bounds=limits_dict['energy_bounds'],**kwargs)
        
    ax.set_ylim(limits_dict['ylim'])
    ax.set_xlim(limits_dict['energy_bounds'])
    
    ax.legend(**legend_dict)
    
    plt.xlabel(axis_dict['label'][0])   
    plt.ylabel(axis_dict['label'][1])
    
    file_name = f"flux_points_{utl.create_roi_name(region_of_interest['name'], region_of_interest['radius_roi'], region_of_interest['e_ref_min'])}"
    
    path_file =  utl.get_path_flux_points() / f'{file_name}{cfg.format_png}'    
    plt.savefig(path_file, bbox_inches='tight')
#     plt.grid(which="both")
    plt.show()
    
    return

