from pydantic import BaseSettings

# Used only to pass path to the file with specifications
# of plot styles
class Settings(BaseSettings):

    PLOT_STYLES_PATH: str = "./plot_styles"
    PLOT_STYLES_FILE: str = "my_plot_style_2.txt"
