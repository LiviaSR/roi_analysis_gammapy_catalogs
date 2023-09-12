#!/usr/bin/env python
# coding: utf-8

# In[ ]:


path_fp_HAWC = "/home/born-again/Documents/GitHub/CTA_projects/flux_points_outside_gammapy_catalogs/HAWC"
path_fp_LHAASO = "/home/born-again/Documents/GitHub/CTA_projects/flux_points_outside_gammapy_catalogs/LHASSO_publishNature"


# In[ ]:


dir_config = "config"
dir_plot_style = "plot_style"
dir_utilities = "utilities"
dir_spectral_models = "spectral_models"
dir_hawc_analysis = "hawc_analysis"
dir_hess_analysis = "hess_analysis"
dir_lhaaso_analysis = "lhaaso_analysis"
dir_cta_simulation = "cta_simulation"
dir_gammapy_catalogs = "gammapy_catalogs"


# In[ ]:


color_lhaaso = "red"
marker_lhaaso = "o"

color_cta = "blue"
marker_cta = "s"


# In[1]:


# catalogs_tags = ["gamma-cat", "hgps", "2hwc", "3hwc", "3fgl", "4fgl", "2fhl", "3fhl"]


# In[ ]:


# markers = ['H', 'D', 'd', 'P', 'X','o', 'v', '^', '<', '>',  '8', 's', 'p', '*', 'h']
markers = ['s','o']
linestyles = ['solid','dotted','dashed','dashdot']


# In[1]:


colors = ["aqua",
"fuchsia",
"peru",
"brown",
"chartreuse",
"chocolate",
"coral",
"khaki",
"darkblue",
"cadetblue",
"pink",
"indigo",
"seagreen",
"crimson",
"darkmagenta",
"orange",
"springgreen",
"plum",
"maroon",
"navy",
"olive",
"skyblue"          
"orange",
"orangered",
"orchid",
"pink",
"plum",
"purple",
"red",
"salmon",
"sienna",
"silver",
"tan",
"teal",
"azure",
"beige",
"ivory",
"black",
"tomato",
"turquoise",
"violet",
"aquamarine",
"wheat",
"white",
"cyan",
"blue",
"lightblue",
"yellowgreen"]


# In[ ]:





# In[ ]:


format_csv  = '.csv'
format_fits = '.fits'
format_dat  = '.dat'
format_png  = '.png'
format_pdf  = '.pdf'
format_svg  = '.svg'
format_yaml = '.yaml'


# In[ ]:





# In[ ]:


sed_type_e2dnde = 'e2dnde'
sed_type_dnde   = 'dnde'


# In[ ]:


unit_deg = 'deg' # Degrees units


# In[ ]:


frame_icrc = "icrs" # International Celestial Reference System (ICRS)


# In[1]:


# dir_analysis = "analysis"
dir_analysis = "analysis"


# In[ ]:


dir_flux_points_tables = "flux_points_tables"


# In[1]:


# parent directories names
dir_flux_points = "flux_points"

dir_tables = "tables"
dir_figures = "figures"
dir_datasets  = "datasets"
dir_models  = "models"

dir_catalogs_roi = "catalogs_roi"

dir_SED = "SED_models"
dir_SED_from_catalogs = "SED_from_catalogs"


# In[ ]:





# In[ ]:


irf_z = [20,40,60]
irf_h = [0.5, 5, 50]
irf_loc = [("cta_north", "North"),("cta_south", "South")]

