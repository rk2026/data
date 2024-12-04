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
    """
    Fetch a file from a GitHub URL
    
    Args:
        url (str): URL of the file to fetch
    
    Returns:
        BytesIO object or None
    """
    try:
        response = requests.get(url)
        response.raise_for_status()
        return BytesIO(response.content)
    except Exception as e:
        st.error(f"Error fetching file: {e}")
        return None

def process_dem_zonal_stats(dem_path, vector_path):
    """
    Perform zonal statistics on DEM raster within vector polygons
    
    Returns:
        GeoDataFrame with zonal statistics for intersecting polygons
    """
    # Read the vector layer
    gdf = gpd.read_file(vector_path)
    
    with rasterio.open(dem_path) as dem:
        raster_bounds = shapely.geometry.box(*dem.bounds)
        
        # Reproject vector to match raster CRS if needed
        if gdf.crs != dem.crs:
            gdf = gdf.to_crs(dem.crs)
        
        # Strictly filter polygons that intersect with raster bounds
        intersecting_gdf = gdf[gdf.intersects(raster_bounds)]
        
        # Initialize results list
        zonal_results = []
        
        # Process each intersecting polygon
        for idx, row in intersecting_gdf.iterrows():
            try:
                # Create a mask for the current polygon
                mask_geometry = [row.geometry]
                
                # Clip raster to polygon
                out_image, out_transform = mask(dem, mask_geometry, crop=True)
                
                # Flatten and remove nodata values
                valid_pixels = out_image[out_image != dem.nodata]
                
                if len(valid_pixels) > 0:
                    # Calculate zonal statistics
                    min_val = float(np.min(valid_pixels))
                    max_val = float(np.max(valid_pixels))
                    
                    # Only include if at least one value is > 0
                    if max_val > 0:
                        # Find pixel locations for min and max
                        min_pixel_idx = np.unravel_index(np.argmin(out_image), out_image.shape)
                        max_pixel_idx = np.unravel_index(np.argmax(out_image), out_image.shape)
                        
                        # Convert pixel indices to geospatial coordinates
                        min_lon, min_lat = rasterio.transform.xy(out_transform, min_pixel_idx[0], min_pixel_idx[1])
                        max_lon, max_lat = rasterio.transform.xy(out_transform, max_pixel_idx[0], max_pixel_idx[1])
                        
                        # Create Point geometries for min and max locations
                        min_point = Point(min_lon, min_lat)
                        max_point = Point(max_lon, max_lat)
                        
                        # Find intersecting polygons for min and max points
                        min_intersecting_polys = intersecting_gdf[intersecting_gdf.contains(min_point)]
                        max_intersecting_polys = intersecting_gdf[intersecting_gdf.contains(max_point)]
                        
                        result_row = {
                            'min_elevation': min_val,
                            'max_elevation': max_val,
                            'min_lon': min_lon,
                            'min_lat': min_lat,
                            'max_lon': max_lon,
                            'max_lat': max_lat,
                            'geometry': row.geometry,
                            'min_point_geometry': min_point,
                            'max_point_geometry': max_point,
                            'min_intersecting_polys': min_intersecting_polys,
                            'max_intersecting_polys': max_intersecting_polys
                        }
                        
                        # Add attributes from the original polygon
                        additional_attrs = ['DISTRICT', 'GaPa_NaPa', 'Type_GN', 'NEW_WARD_N']
                        for attr in additional_attrs:
                            if attr in row.index:
                                result_row[attr] = row[attr]
                        
                        zonal_results.append(result_row)
            
            except Exception as e:
                st.warning(f"Could not process polygon: {e}")
        
        # Create GeoDataFrame from results
        results_gdf = gpd.GeoDataFrame(zonal_results, crs=gdf.crs)
        return results_gdf, intersecting_gdf

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
