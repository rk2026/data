def create_3d_visualization(dem_path, intersecting_gdf, zonal_results):
    """
    Create enhanced 3D visualization combining DEM and vector layers
    
    Returns:
        Plotly Figure
    """
    # Read the DEM raster
    with rasterio.open(dem_path) as dem:
        # Read raster data
        dem_array = dem.read(1)
        
        # Get raster bounds and transform
        bounds = dem.bounds
        
        # Create x and y coordinates
        x = np.linspace(bounds.left, bounds.right, dem_array.shape[1])
        y = np.linspace(bounds.bottom, bounds.top, dem_array.shape[0])
        
        # Calculate z-scale factor (typically < 1 to avoid vertical exaggeration)
        elevation_range = dem_array.max() - dem_array.min()
        lat_range = bounds.top - bounds.bottom
        lon_range = bounds.right - bounds.left
        z_scale = min(lat_range, lon_range) / elevation_range * 0.1
        
        # Create 3D surface plot of DEM with adjusted z-scale
        surface_trace = go.Surface(
            z=dem_array * z_scale, 
            x=x, 
            y=y, 
            colorscale='Viridis', 
            showscale=False,  # Remove color scale legend
            name='Terrain Elevation',
            opacity=0.7
        )
        
        # Prepare data for traces
        traces = [surface_trace]
        
        # Add vector layer boundaries
        for _, polygon in intersecting_gdf.iterrows():
            # Extract exterior coordinates of the polygon
            if polygon.geometry.type == 'Polygon':
                coords = list(polygon.geometry.exterior.coords)
            elif polygon.geometry.type == 'MultiPolygon':
                # For multipolygon, use the first polygon's exterior
                coords = list(list(polygon.geometry.geoms)[0].exterior.coords)
            
            # Get z values for the polygon boundary
            boundary_z = np.interp(
                [coord[2] if len(coord) > 2 else dem_array.mean() for coord in coords], 
                [dem_array.min(), dem_array.max()], 
                [dem_array.min() * z_scale, dem_array.max() * z_scale]
            )
            
            # Create polygon boundary trace
            boundary_trace = go.Scatter3d(
                x=[coord[0] for coord in coords],
                y=[coord[1] for coord in coords],
                z=boundary_z,
                mode='lines',
                line=dict(color='red', width=2),
                showlegend=False
            )
            traces.append(boundary_trace)
        
        # Add min and max elevation points with draped positioning
        for _, row in zonal_results.iterrows():
            # Find the closest grid point for min elevation
            min_x_idx = np.argmin(np.abs(x - row['min_lon']))
            min_y_idx = np.argmin(np.abs(y - row['min_lat']))
            min_z = dem_array[min_y_idx, min_x_idx] * z_scale
            
            # Find the closest grid point for max elevation
            max_x_idx = np.argmin(np.abs(x - row['max_lon']))
            max_y_idx = np.argmin(np.abs(y - row['max_lat']))
            max_z = dem_array[max_y_idx, max_x_idx] * z_scale
            
            # Minimum elevation point (draped)
            min_trace = go.Scatter3d(
                x=[row['min_lon']], 
                y=[row['min_lat']], 
                z=[min_z],
                mode='markers',
                marker=dict(
                    size=10,
                    color='blue',
                    opacity=0.8
                ),
                showlegend=False,
                text=f"Min Elev: {row.get('min_elevation', 'N/A')}m<br>"
                     f"Ward: {row.get('NEW_WARD_N', 'N/A')}"
            )
            traces.append(min_trace)
            
            # Maximum elevation point (draped)
            max_trace = go.Scatter3d(
                x=[row['max_lon']], 
                y=[row['max_lat']], 
                z=[max_z],
                mode='markers',
                marker=dict(
                    size=10,
                    color='red',
                    opacity=0.8
                ),
                showlegend=False,
                text=f"Max Elev: {row.get('max_elevation', 'N/A')}m<br>"
                     f"Ward: {row.get('NEW_WARD_N', 'N/A')}"
            )
            traces.append(max_trace)
        
        # Create 3D figure
        fig = go.Figure(data=traces)
        
        # Update layout for full-screen and better view
        fig.update_layout(
            title='3D Terrain Visualization',
            scene=dict(
                xaxis_title='Longitude',
                yaxis_title='Latitude',
                zaxis_title='Elevation (scaled)',
                aspectmode='manual',
                aspectratio=dict(x=1, y=1, z=0.3),  # Flatten z-axis
                camera=dict(
                    eye=dict(x=1.5, y=1.5, z=1)
                )
            ),
            height=900,  # Increased height
            width=1600,  # Significantly increased width
            margin=dict(l=0, r=0, t=30, b=0),  # Reduced margins
            showlegend=False  # Remove legend completely
        )
        
        return fig
