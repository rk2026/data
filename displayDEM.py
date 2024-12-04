import streamlit as st
import plotly.express as px
import geopandas as gpd
import rasterio
import numpy as np
import plotly.graph_objs as go
from rasterio.mask import mask
import requests
from io import BytesIO
import shapely
from shapely.geometry import Point

# Optional imports with error handling
try:
    import osmnx as ox
    OSM_AVAILABLE = True
except ImportError:
    st.warning("OpenStreetMap (osmnx) libraries not installed. OSM features will be disabled.")
    OSM_AVAILABLE = False

DEM_FILES = [
    'DHADING_Netrawati.tif',
    'DHADING_Khaniyabash.tif', 
    'DHADING_Jwalamukhi.tif', 
    'DHADING_Galchi.tif', 
    'DHADING_Gajuri.tif', 
    'DHADING_Benighat_Rorang.tif', 
    'DHADING_Nilkantha.tif',
    'DHADING_Gangajamuna.tif',
    'DHADING_Dhunibesi.tif',
    'DHADING_Siddhalek.tif',
]

def fetch_github_file(url):
    try:
        response = requests.get(url)
        response.raise_for_status()
        return BytesIO(response.content)
    except Exception as e:
        st.error(f"Error fetching file: {e}")
        return None

# ... [previous functions remain the same] ...

def create_2d_map(intersecting_gdf, zonal_results):
    """
    Create a 2D map visualization with vector layers and elevation points
    
    Returns:
        Plotly Figure
    """
    # Calculate center of the map
    center_lon = intersecting_gdf.geometry.centroid.x.mean()
    center_lat = intersecting_gdf.geometry.centroid.y.mean()
    
    # Create base map
    fig = go.Figure()
    
    # Add Ward Boundaries with labels
    for _, polygon in intersecting_gdf.iterrows():
        # Determine polygon coordinates
        if polygon.geometry.type == 'Polygon':
            coords = list(polygon.geometry.exterior.coords)
            lons = [coord[0] for coord in coords]
            lats = [coord[1] for coord in coords]
        elif polygon.geometry.type == 'MultiPolygon':
            # For multipolygon, use the first polygon's exterior
            coords = list(list(polygon.geometry.geoms)[0].exterior.coords)
            lons = [coord[0] for coord in coords]
            lats = [coord[1] for coord in coords]
        
        # Add polygon boundary
        fig.add_trace(
            go.Scattermapbox(
                mode='lines',
                lon=lons,
                lat=lats,
                line=dict(color='red', width=2),
                showlegend=False
            )
        )
        
        # Add polygon centroid label
        centroid = polygon.geometry.centroid
        label_text = f"{polygon.get('GaPa_NaPa', 'N/A')}<br>Ward: {polygon.get('NEW_WARD_N', 'N/A')}"
        
        fig.add_trace(
            go.Scattermapbox(
                mode='markers+text',
                lon=[centroid.x],
                lat=[centroid.y],
                marker=dict(
                    size=8,
                    color='green',
                    opacity=0.5
                ),
                text=[label_text],
                textposition='bottom center',
                hoverinfo='text',
                showlegend=False
            )
        )
    
    # Add Minimum and Maximum Elevation Points
    for _, row in zonal_results.iterrows():
        # Minimum Elevation Point
        # Get the intersecting polygon details for min point
        min_poly_details = row.get('min_intersecting_polys', [])
        min_poly_info = min_poly_details.iloc[0] if len(min_poly_details) > 0 else None
        
        fig.add_trace(
            go.Scattermapbox(
                mode='markers',
                lon=[row['min_lon']],
                lat=[row['min_lat']],
                marker=dict(
                    size=10,
                    color='blue',
                    opacity=0.7
                ),
                text=f"Minimum Elevation: {row.get('min_elevation', 'N/A')}m<br>"
                     f"Ward: {min_poly_info.get('NEW_WARD_N', 'N/A') if min_poly_info is not None else 'N/A'}<br>"
                     f"District: {min_poly_info.get('DISTRICT', 'N/A') if min_poly_info is not None else 'N/A'}",
                hoverinfo='text',
                showlegend=False
            )
        )
        
        # Maximum Elevation Point
        # Get the intersecting polygon details for max point
        max_poly_details = row.get('max_intersecting_polys', [])
        max_poly_info = max_poly_details.iloc[0] if len(max_poly_details) > 0 else None
        
        fig.add_trace(
            go.Scattermapbox(
                mode='markers',
                lon=[row['max_lon']],
                lat=[row['max_lat']],
                marker=dict(
                    size=10,
                    color='red',
                    opacity=0.7
                ),
                text=f"Maximum Elevation: {row.get('max_elevation', 'N/A')}m<br>"
                     f"Ward: {max_poly_info.get('NEW_WARD_N', 'N/A') if max_poly_info is not None else 'N/A'}<br>"
                     f"District: {max_poly_info.get('DISTRICT', 'N/A') if max_poly_info is not None else 'N/A'}",
                hoverinfo='text',
                showlegend=False
            )
        )
    
    # Update layout for OpenStreetMap style with map type toggle
    fig.update_layout(
        mapbox_style="open-street-map",
        mapbox=dict(
            center=dict(
                lat=center_lat,
                lon=center_lon
            ),
            zoom=9
        ),
        height=800,
        width=1200,
        title='Maps of Minimum and Maximum Location Elevations',
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                x=0.57,
                y=1.2,
                showactive=True,
                buttons=list([
                    dict(label="Street Map",
                         method="relayout",
                         args=[{"mapbox.style": "open-street-map"}]),
                    dict(label="Satellite",
                         method="relayout", 
                         args=[{"mapbox.style": "satellite-streets"}]),
                    dict(label="Terrain",
                         method="relayout", 
                         args=[{"mapbox.style": "carto-positron"}])
                ]),
            )
        ],
        margin={"r":0,"t":50,"l":0,"b":0}
    )
    
    return fig

# ... [rest of the code remains the same] ...

def main():
    st.title("DEM Analysis and Visualization")
    
    # Dropdown for selecting DEM file
    selected_dem = st.selectbox("Select DEM File", DEM_FILES)
    
    # Construct GitHub URL for the selected DEM
    base_url = "https://github.com/rk2026/data/raw/main/"
    dem_url = base_url + selected_dem
    
    # Fixed Vector Layer URL
    vector_url = base_url + "Bagmati_ward.gpkg"
    
    if st.button("Process Data"):
        with st.spinner(f"Processing {selected_dem}..."):
            # Fetch files
            dem_file = fetch_github_file(dem_url)
            vector_file = fetch_github_file(vector_url)
            
            if dem_file and vector_file:
                # Process zonal statistics
                zonal_results, intersecting_gdf = process_dem_zonal_stats(dem_file, vector_file)
                
                # Display results table
                st.subheader(f"Zonal Statistics Results for {selected_dem}")
                # Select columns to display
                display_columns = [
                    'DISTRICT', 'GaPa_NaPa', 'Type_GN', 'NEW_WARD_N', 
                    'min_elevation', 'max_elevation', 
                    'min_lon', 'min_lat', 'max_lon', 'max_lat'
                ]
                
                # Filter columns that exist in the results
                available_columns = [col for col in display_columns if col in zonal_results.columns]
                
                st.dataframe(zonal_results[available_columns])
                
                # Verify and display number of intersecting polygons
                st.write(f"Number of intersecting polygons: {len(intersecting_gdf)}")
                st.write(f"Number of zonal results: {len(zonal_results)}")
                
                # Additional details about point intersections
                st.write("Minimum Point Intersections:")
                for idx, row in zonal_results.iterrows():
                    min_polys = row.get('min_intersecting_polys', [])
                    max_polys = row.get('max_intersecting_polys', [])
                    
                    st.write(f"Result {idx}:")
                    st.write(f"  Minimum Point Intersecting Polygons: {len(min_polys)}")
                    st.write(f"  Maximum Point Intersecting Polygons: {len(max_polys)}")
                
                # Create and display 2D map
                fig_2d = create_2d_map(intersecting_gdf, zonal_results)
                st.plotly_chart(fig_2d, use_container_width=True)

if __name__ == "__main__":
    main()
