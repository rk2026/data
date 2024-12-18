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
    Calculate zonal statistics for each polygon in the vector layer
    """
    try:
        # Read vector layer
        vector_gdf = gpd.read_file(vector_path)
        
        # Read the raster and vector data
        with rasterio.open(dem_path) as dem_src:
            # Ensure consistent CRS
            vector_gdf = vector_gdf.to_crs(dem_src.crs)
            
            # Get raster bounds
            raster_bounds = box(*dem_src.bounds)
            
            # Initialize lists to store results
            results = []
            
            # Iterate through each polygon
            for index, row in vector_gdf.iterrows():
                try:
                    # Check if polygon intersects with raster bounds
                    if not row.geometry.intersects(raster_bounds):
                        st.warning(f"Polygon {index} does not intersect with the raster.")
                        continue
                    
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
                        
                        # Ensure we have valid pixel coordinates
                        if min_pixel[0].size > 0 and min_pixel[1].size > 0 and \
                           max_pixel[0].size > 0 and max_pixel[1].size > 0:
                            # Convert pixel coordinates to geographic coordinates
                            min_lon, min_lat = rasterio.transform.xy(
                                out_transform, min_pixel[1][0], min_pixel[0][0]
                            )
                            max_lon, max_lat = rasterio.transform.xy(
                                out_transform, max_pixel[1][0], max_pixel[0][0]
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
                                'max_lat': max_lat,
                                'geometry': row.geometry
                            })
                        else:
                            st.warning(f"Could not find valid pixel coordinates for polygon {index}")
                
                except Exception as e:
                    st.warning(f"Could not process polygon {index}: {str(e)}")
            
            # Convert results to DataFrame
            results_df = pd.DataFrame(results)
            return results_df, vector_gdf
    
    except Exception as e:
        st.error(f"Error in zonal statistics calculation: {e}")
        return pd.DataFrame(), gpd.GeoDataFrame()

def create_interactive_map(dem_results, vector_gdf):
    """
    Create an interactive map with administrative boundaries and points
    """
    try:
        # Print column names for debugging
        print("Columns in dem_results:", list(dem_results.columns))
        
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
            # Check if required columns exist
            required_columns = ['max_lat', 'max_lon', 'min_lat', 'min_lon', 'min_elevation', 'max_elevation', 'NEW_WARD_N']
            if not all(col in dem_results.columns for col in required_columns):
                st.warning(f"Missing required columns for mapping. Available columns: {list(dem_results.columns)}")
                continue
            
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
