import streamlit as st
import geopandas as gpd
import requests
from io import BytesIO
from streamlit_folium import st_folium
import folium

# Set page title
st.title("Display GeoPackage from GitHub")

# GitHub file URL (replace with your raw URL)
github_url = "https://raw.githubusercontent.com/rk2026/data/main/topo_grid.gpkg"

# Function to fetch the GeoPackage
@st.cache_data
def fetch_geopackage(url):
    response = requests.get(url)
    if response.status_code == 200:
        return BytesIO(response.content)
    else:
        st.error("Failed to fetch the file. Check the URL.")
        return None

# Fetch the file
file_content = fetch_geopackage(github_url)

if file_content:
    # Load the GeoPackage into a GeoDataFrame
    try:
        gdf = gpd.read_file(file_content)
        st.success("GeoPackage loaded successfully!")
        
        # Display GeoDataFrame info
        st.subheader("GeoPackage Info")
        st.write(gdf.head())
        
        # Display map using Folium
        st.subheader("Interactive Map")
        centroid = gdf.geometry.unary_union.centroid
        m = folium.Map(location=[centroid.y, centroid.x], zoom_start=10)
        folium.GeoJson(gdf).add_to(m)
        st_folium(m, width=800, height=500)
    except Exception as e:
        st.error(f"Error loading GeoPackage: {e}")
