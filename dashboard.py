import os
import sys
import streamlit as st

from matplotlib import pyplot as plt, patches

from gammapy.datasets import Datasets
from gammapy.maps import RegionGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    Models, 
    SkyModel, 
    PowerLawSpectralModel,
    LogParabolaSpectralModel,
)
from gammapy.utils.check import check_tutorials_setup
from gammapy.visualization.utils import plot_contour_line

from astropy import units as u
from astropy.coordinates import SkyCoord 
from regions import CircleSkyRegion, PointSkyRegion

import my_modules.config.cfg as cfg
import my_modules.plot_style.plotter as plotter
import my_modules.utilities.utilities as utl
from my_modules.spectral_models import SkyModelFactory
import my_modules.gammapy_catalogs.gammapy_catalogs as gammapy_cat

from settings import Settings

settings = Settings()





# Sort systems by name
def sort_by_system_name(systems):
    sorted_systems = sorted(systems, key=lambda s: s[1:])
    return sorted_systems


# Function to generate figures
def list_images():
    images = []
    for system_name in os.listdir("./roi_analysis_output"):
        if system_name == ".DS_Store": continue
        for image in os.listdir(f"./roi_analysis_output/{system_name}"):
            if image[-3:] == "png":
                images.append(f"{system_name}/{image}")
    return sort_by_system_name(images)

# Streamlit App
def main():

# Plot the sky map of the region of interest
    st.markdown("## Map Region of Interest")

    radius_roi = st.slider('RoI size', 0.5, 2.0, 0.5, 0.5)

    path_my_plot_style = f"{settings.PLOT_STYLES_PATH}/{settings.PLOT_STYLES_FILE}" 
    plt.style.use(path_my_plot_style)

    filename = "SpiderSystems.csv"

    spider_systems = utl.read_systems_file(filename)

    selected_spider_system = st.selectbox("Select System", [spider_system.source_name for spider_system in spider_systems])

    system_data = list(filter(lambda e: e.source_name == selected_spider_system, spider_systems))[0]
    print(system_data)

    source_info = utl.set_source_info(system_data)

    region_of_interest = utl.create_region_of_interest(
    source_info=source_info, 
    radius_roi=radius_roi, 
    e_ref_min=None, 
    e_ref_max=None,
    )
    
    sources_gammapy, _, _ = gammapy_cat.get_datasets_flux_points_gammapy(region_of_interest)

    # Loop to set markers and color of each source
    while len(cfg.markers) < len(sources_gammapy) +1:
        cfg.markers.extend(cfg.markers)
    while len(cfg.colors) < len(sources_gammapy) +1:
        cfg.colors.extend(cfg.colors)

    # Set a circle (region of interest) around the source

    center = region_of_interest["position"]
    angle = radius_roi * u.deg

    circle = RegionGeom(CircleSkyRegion(center, angle))
    fig = plt.figure(figsize=(6,6))
    ax = circle.plot_region()

    point = RegionGeom(PointSkyRegion(center=center))
    ax = point.plot_region(
        ax=ax,
        facecolor="black",
        edgecolor="black",
        kwargs_point={
            "color": "black",
            "fillstyle": "full",
            "marker": "X",
            "markersize": 15,
        },
    )
    #print(point.axes)

    for i, source in enumerate(sources_gammapy):
        label = plotter.set_label_datasets(source.name)
        point_cat = RegionGeom(PointSkyRegion(center=source.position))
        point_cat.plot_region(
            ax=ax,
            label=label,
            facecolor=cfg.colors[i],
            edgecolor=cfg.colors[i],
            kwargs_point={
                "color": cfg.colors[i],
                "fillstyle": "full",
                "marker": cfg.markers[i],
                "markersize": 10,
            },
        )


    _, center, _ = st.columns([1,10, 1])
    with center:
        st.pyplot(fig, clear_figure=True)

#######################################################################################

    st.title("Flux Points Spider Systems")

    # Dropdown menu to select the figure type
    figure = st.selectbox("Select System", list_images())

    # Button to generate the selected figure
    if st.button("Generate Figure"):
        path = f"./roi_analysis_output/{figure}"

        # Display the figure
        st.image(path)
    
    # Add a section to display multiple figures
    st.markdown("## Display Multiple Figures")

    # Get a list of all available systems
    all_systems = list_images()

        # Select multiple systems
    selected_systems = st.multiselect("Select Systems to Display", all_systems)

        # Button to generate the selected figures
    if st.button("Generate Selected Figures"):
        # Display the selected figures side by side
        num_selected = len(selected_systems)
#    rows = len(selected_systems) // 2
#    columns = st.columns(2)
#    for system, column in zip(selected_systems, columns):
#        with column:
#            path = f"./roi_analysis_output/{system}"
#            st.image(path, caption=system)
        rows = (num_selected + 1) // 2
        columns = 2

        for i in range(rows):
            col1, col2 = st.columns(2)
            for j in range(columns):
                index = i * columns + j
                if index < num_selected:
                    image_path = os.path.join("roi_analysis_output", selected_systems[index])
                    if j == 0:
                        col1.image(image_path, caption=selected_systems[index], use_column_width=True)
                    else:
                        col2.image(image_path, caption=selected_systems[index], use_column_width=True)


if __name__ == "__main__":
    main()

