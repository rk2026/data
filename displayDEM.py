import streamlit as st
import geopandas as gpd
import rasterio
from rasterio.mask import mask
import numpy as np
import pandas as pd
import folium
from folium import GeoJson
from streamlit_folium import folium_static
from shapely.geometry import box

def load_dem_files():
    """
    Load predefined DEM files from GitHub or local sources
    """
    dem_files = {
        'DHADING_Thakre.tif': 'https://raw.githubusercontent.com/rk2026/data/main/DHADING_Thakre.tif',
        'DHADING_Gajuri.tif': 'https://raw.githubusercontent.com/rk2026/data/main/DHADING_Gajuri.tif',
    }
    vector_path = 'https://raw.githubusercontent.com/rk2026/data/main/Bagmati_ward.gpkg'
    return dem_files, vector_path

def calculate_zonal_statistics(dem_path, vector_path):
    """
    Calculate zonal statistics with robust spatial intersection handling
    """
    try:
        # Read vector layer
        vector_gdf = gpd.read_file(vector_path)
        
        # Read the raster and vector data
        with rasterio.open(dem_path) as dem_src:
            # Ensure consistent CRS
            vector_gdf = vector_gdf.to_crs(dem_src.crs)
            
            # Get raster bounds and create spatial index for faster intersection
            raster_bounds = box(*dem_src.bounds)
            
            # Initialize lists to store results
            results = []
            
            # Iterate through each polygon
            for index, row in vector_gdf.iterrows():
                try:
                    # Precise intersection check with actual geometry
                    intersection = row.geometry.intersection(raster_bounds)
                    
                    # Skip if no meaningful intersection
                    if intersection.is_empty or intersection.area == 0:
                        st.warning(f"Polygon {index} does not meaningfully intersect with the raster.")
                        continue
                    
                    # Clip polygon to raster bounds
                    clipped_geom = row.geometry.intersection(raster_bounds)
                    
                    # Create a mask for the clipped geometry
                    out_image, out_transform = mask(
                        dem_src, 
                        [clipped_geom.__geo_interface__], 
                        crop=True, 
                        nodata=np.nan
                    )
                    
                    # Rest of your existing statistics calculation...
                    masked_data = out_image[0]
                    valid_data = masked_data[~np.isnan(masked_data)]
                    
                    if len(valid_data) > 0:
                        min_val = np.min(valid_data)
                        max_val = np.max(valid_data)
                        
                        # Existing coordinate and result processing...
                        results.append({
                            'DISTRICT': row.get('DISTRICT', 'N/A'),
                            'GaPa_NaPa': row.get('GaPa_NaPa', 'N/A'),
                            'NEW_WARD_N': row.get('NEW_WARD_N', 'N/A'),
                            'Type_GN': row.get('Type_GN', 'N/A'),
                            'min_elevation': min_val,
                            'max_elevation': max_val,
                            # ... other existing fields
                        })
                
                except Exception as e:
                    st.warning(f"Error processing polygon {index}: {str(e)}")
            
            results_df = pd.DataFrame(results)
            return results_df, vector_gdf
    
    except Exception as e:
        st.error(f"Comprehensive error in zonal statistics: {e}")
        return pd.DataFrame(), gpd.GeoDataFrame()

def create_interactive_map(dem_results, vector_gdf):
    """
    Create an interactive map with administrative boundaries and points
    """
    try:
        # Compute map center
        center_lat = vector_gdf.geometry.centroid.y.mean()
        center_lon = vector_gdf.geometry.centroid.x.mean()
        
        # Create map
        m = folium.Map(
            location=[center_lat, center_lon], 
            zoom_start=8,
            tiles='OpenStreetMap'
        )
        
        # Add basemap layers
        folium.TileLayer(
            tiles='OpenStreetMap',
            name='OpenStreetMap',
            attr='OpenStreetMap Contributors'
        ).add_to(m)
        
        folium.TileLayer(
            tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
            attr='Esri',
            name='Satellite',
            overlay=False
        ).add_to(m)
        
        # Add ward boundary layer
        ward_style = {
            'fillColor': 'blue',
            'color': 'black',
            'weight': 2,
            'fillOpacity': 0.1
        }
        
        def ward_popup(feature):
            props = feature['properties']
            return folium.Popup(f"""
                District: {props.get('DISTRICT', 'N/A')}
                Ward: {props.get('NEW_WARD_N', 'N/A')}
                Type: {props.get('Type_GN', 'N/A')}
            """)
        
        # Add ward boundaries
        GeoJson(
            vector_gdf.__geo_interface__,
            name='Ward Boundaries',
            style_function=lambda x: ward_style,
            popup=ward_popup
        ).add_to(m)
        
        # Add minimum and maximum elevation points
        for idx, row in dem_results.iterrows():
            # Minimum height point
            folium.CircleMarker(
                location=[row['min_lat'], row['min_lon']],
                radius=6,
                popup=f"Minimum Elevation: {row['min_elevation']:.2f}m\nWard: {row['NEW_WARD_N']}",
                color='green',
                fill=True,
                fillColor='green',
                fillOpacity=0.7
            ).add_to(m)
            
            # Maximum height point
            folium.CircleMarker(
                location=[row['max_lat'], row['max_lon']],
                radius=6,
                popup=f"Maximum Elevation: {row['max_elevation']:.2f}m\nWard: {row['NEW_WARD_N']}",
                color='red',
                fill=True,
                fillColor='red',
                fillOpacity=0.7
            ).add_to(m)
        
        # Add layer control
        folium.LayerControl().add_to(m)
        
        return m
    except Exception as e:
        st.error(f"Error creating interactive map: {e}")
        return None

def main():
    st.set_page_config(page_title="DEM Zonal Statistics Analysis", layout="wide")
    
    st.title('DEM Zonal Statistics Analysis')
    
    # Sidebar for file selection
    st.sidebar.header('Data Selection')
    
    # DEM and Vector File Selection
    dem_files, vector_path = load_dem_files()
    
    # Add file upload option
    uploaded_dem = st.sidebar.file_uploader(
        "Upload DEM File", 
        type=['.tif', '.tiff']
    )
    
    # Combine predefined and uploaded files
    if uploaded_dem:
        dem_files['Uploaded DEM'] = uploaded_dem
    
    selected_dem = st.sidebar.selectbox(
        'Select DEM File', 
        list(dem_files.keys())
    )
    
    # Process data
    if st.sidebar.button('Analyze DEM'):
        with st.spinner('Processing DEM and calculating zonal statistics...'):
            try:
                # Determine the path to use
                if uploaded_dem and selected_dem == 'Uploaded DEM':
                    dem_path = uploaded_dem
                else:
                    dem_path = dem_files[selected_dem]
                
                # Calculate zonal statistics
                dem_results, vector_gdf = calculate_zonal_statistics(dem_path, vector_path)
                
                if not dem_results.empty:
                    # Create two columns for display
                    col1, col2 = st.columns(2)
                    
                    # Display results table in the first column
                    with col1:
                        st.subheader("Zonal Statistics Results")
                        st.dataframe(dem_results[['DISTRICT', 'GaPa_NaPa', 'NEW_WARD_N', 'min_elevation', 'max_elevation']])
                    
                    # Create interactive map in the second column
                    with col2:
                        st.subheader("Interactive Map")
                        map_obj = create_interactive_map(dem_results, vector_gdf)
                        if map_obj:
                            folium_static(map_obj)
                        else:
                            st.warning("Could not create interactive map")
                else:
                    st.warning("No valid results found")
                
            except Exception as e:
                st.error(f"An error occurred: {str(e)}")

if __name__ == '__main__':
    main()
