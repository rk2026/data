import streamlit as st
import rasterio
from rasterio.plot import show
import folium
from streamlit_folium import st_folium
from folium.raster_layers import ImageOverlay
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize
import requests
from io import BytesIO

# Set Streamlit app title
st.title("Display DEM File from GitHub")

# GitHub file URL (update this with your file's raw URL)
github_url = "https://raw.githubusercontent.com/rk2026/data/main/DHADING_Thakre.tif"

# Load DEM file from GitHub
try:
    st.info("Downloading the DEM file from GitHub...")
    response = requests.get(github_url)
    response.raise_for_status()
    dem_file = BytesIO(response.content)
    
    # Read the DEM file
    with rasterio.open(dem_file) as src:
        dem_data = src.read(1)  # Read the first band
        dem_bounds = src.bounds
        dem_transform = src.transform
    
    # Normalize DEM data for visualization
    min_val, max_val = np.min(dem_data), np.max(dem_data)
    norm = Normalize(vmin=min_val, vmax=max_val)
    cmap = cm.terrain
    colormap = cmap(norm(dem_data))
    
    # Convert colormap to RGBA image
    rgba_image = (colormap * 255).astype(np.uint8)
    
#     folium_map = folium.Map(
#     location=[(dem_bounds.top + dem_bounds.bottom) / 2, (dem_bounds.left + dem_bounds.right) / 2],
#     zoom_start=12,
#     tiles="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
#     attr='Google Satellite',
# )

    
    # Create Folium Map
   folium_map = folium.Map(
   location=[(dem_bounds.top + dem_bounds.bottom) / 2, (dem_bounds.left + dem_bounds.right) / 2],
   zoom_start=12,
   tiles="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
   attr='Google Satellite',
)
    
    # Add Image Overlay to Folium Map
    image_overlay = ImageOverlay(
        rgba_image,
        bounds=[[dem_bounds.bottom, dem_bounds.left], [dem_bounds.top, dem_bounds.right]],
        opacity=0.8,
    )
    image_overlay.add_to(folium_map)
    
    # Display the Folium Map in Streamlit
    st_folium(folium_map, width=700, height=500)
except requests.exceptions.RequestException as e:
    st.error(f"Failed to download the file: {e}")
