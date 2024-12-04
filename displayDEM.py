import streamlit as st
import geopandas as gpd
import rasterio
from rasterio.mask import mask
import numpy as np
import pandas as pd
import folium
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

def filter_intersecting_polygons(vector_gdf, raster_path):
    """
    Filter vector polygons that intersect with the raster extent
    """
    try:
        # Open the raster to get its bounds and CRS
        with rasterio.open(raster_path) as raster_src:
            raster_crs = raster_src.crs
            raster_bounds = raster_src.bounds
        
        # Reproject vector layer to match the raster CRS
        vector_gdf = vector_gdf.to_crs(raster_crs)
        
        # Create a bounding box from the raster extent
        raster_bbox = box(raster_bounds.left, raster_bounds.bottom, 
                          raster_bounds.right, raster_bounds.top)
        
        # Filter polygons that intersect with the raster extent
        intersecting_gdf = vector_gdf[vector_gdf.geometry.intersects(raster_bbox)]
        
        return intersecting_gdf
    except Exception as e:
        st.error(f"Error filtering intersecting polygons: {e}")
        return gpd.GeoDataFrame()

def calculate_zonal_statistics(dem_path, vector_path):
    """
    Calculate zonal statistics for each polygon in the vector layer
    """
    try:
        # Read vector layer
        vector_gdf = gpd.read_file(vector_path)
        
        # Filter polygons that intersect with raster extent
        vector_gdf = filter_intersecting_polygons(vector_gdf, dem_path)
        
        # Read the raster and vector data
        with rasterio.open(dem_path) as dem_src:
            # Initialize lists to store results
            results = []
            
            # Iterate through each polygon
            for index, row in vector_gdf.iterrows():
                try:
                    # Create a mask for the current polygon
                    out_image, out_transform = mask(
                        dem_src, 
                        [row.geometry.__geo_interface__], 
                        crop=True, 
                        nodata=np.nan
                    )
                    
                    # Calculate statistics
                    masked_data = out_image[0]
                    valid_data = masked_data[~np.isnan(masked_data)]
                    
                    if len(valid_data) > 0:
                        min_val = np.min(valid_data)
                        max_val = np.max(valid_data)
                        
                        # Find pixel coordinates for min and max
                        min_pixel = np.where(masked_data == min_val)
                        max_pixel = np.where(masked_data == max_val)
                        
                        # Convert pixel coordinates to geographic coordinates
                        min_lon, min_lat = rasterio.transform.xy(
                            out_transform, min_pixel[0][0], min_pixel[1][0]
                        )
                        max_lon, max_lat = rasterio.transform.xy(
                            out_transform, max_pixel[0][0], max_pixel[1][0]
                        )
                        
                        results.append({
                            'DISTRICT': row.get('DISTRICT', 'N/A'),
                            'GaPa_NaPa': row.get('GaPa_NaPa', 'N/A'),
                            'NEW_WARD_N': row.get('NEW_WARD_N', 'N/A'),
                            'Type_GN': row.get('Type_GN', 'N/A'),
                            'min_elevation': min_val,
                            'max_elevation': max_val,
                            'min_lon': min_lon,
                            'min_lat': min_lat,
                            'max_lon': max_lon,
                            'max_lat': max_lat
                        })
                except Exception as e:
                    st.warning(f"Could not process polygon {index}: {str(e)}")
            
            # Convert results to DataFrame
            results_df = pd.DataFrame(results)
            return results_df
    except Exception as e:
        st.error(f"Error in zonal statistics calculation: {e}")
        return pd.DataFrame()

def create_interactive_map(dem_results, vector_path):
    """
    Create an interactive map with administrative boundaries and points
    """
    try:
        # Read vector layer
        vector_gdf = gpd.read_file(vector_path)
        
        # Create a map centered on the mean coordinates of the vector layer
        center_lat = vector_gdf.geometry.centroid.y.mean()
        center_lon = vector_gdf.geometry.centroid.x.mean()
        
        # Create map with explicit attribution
        m = folium.Map(
            location=[center_lat, center_lon], 
            zoom_start=10,
            tiles='OpenStreetMap',
            attr='&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
        )
        
        # Add additional base layers with explicit attribution
        folium.TileLayer(
            tiles='https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png',
            attr='&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
            name='OpenStreetMap'
        ).add_to(m)
        
        folium.TileLayer(
            tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
            attr='Tiles &copy; Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community',
            name='Satellite'
        ).add_to(m)
        
        # Add layer control
        folium.LayerControl().add_to(m)
        
        # Style function for vector layer
        def style_function(feature):
            return {
                'fillColor': 'blue',
                'color': 'black',
                'weight': 1,
                'fillOpacity': 0.1
            }
        
        # Popup function for vector layer
        def popup_function(feature):
            props = feature['properties']
            return folium.Popup(f"""
                District: {props.get('DISTRICT', 'N/A')}
                Ward: {props.get('NEW_WARD_N', 'N/A')}
                Type: {props.get('Type_GN', 'N/A')}
            """)
        
        # Add vector layer
        folium.GeoJson(
            vector_gdf, 
            style_function=style_function,
            popup=popup_function
        ).add_to(m)
        
        # Add minimum and maximum height points
        for idx, row in dem_results.iterrows():
            # Minimum height point
            folium.CircleMarker(
                location=[row['min_lat'], row['min_lon']],
                radius=5,
                popup=f"Minimum Elevation: {row['min_elevation']:.2f}m",
                color='green',
                fill=True,
                fillColor='green'
            ).add_to(m)
            
            # Maximum height point
            folium.CircleMarker(
                location=[row['max_lat'], row['max_lon']],
                radius=5,
                popup=f"Maximum Elevation: {row['max_elevation']:.2f}m",
                color='red',
                fill=True,
                fillColor='red'
            ).add_to(m)
        
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
                dem_results = calculate_zonal_statistics(dem_path, vector_path)
                
                if dem_results is not None and not dem_results.empty:
                    # Create two columns for display
                    col1, col2 = st.columns(2)
                    
                    # Display results table in the first column
                    with col1:
                        st.subheader("Zonal Statistics Results")
                        st.dataframe(dem_results)
                    
                    # Create interactive map in the second column
                    with col2:
                        st.subheader("Interactive Map")
                        map_obj = create_interactive_map(dem_results, vector_path)
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
