import geopandas as gpd
import pandas as pd
import os
import numpy as np

def diagnose():
    shp_path = "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/training_punten_GEE_2.shp"
    downloads_dir = "downloads"
    
    # 1. Load SHP
    print("Loading SHP...")
    try:
        gdf = gpd.read_file(shp_path)
        print(f"Columns: {gdf.columns.tolist()}")
        
        # Identify columns
        lat_col = next((c for c in gdf.columns if 'lat' in c.lower()), None)
        lon_col = next((c for c in gdf.columns if 'lon' in c.lower()), None)
        soil_col = next((c for c in gdf.columns if 'soil' in c.lower() or 'pure' in c.lower()), None)
        
        print(f"Using: {lat_col}, {lon_col}, {soil_col}")
        
        pure_points = []
        if soil_col:
            # Check what "pure" means
            print(f"Values in {soil_col}: {gdf[soil_col].unique()}")
            
            for idx, row in gdf.iterrows():
                val = row[soil_col]
                is_pure = False
                try:
                    if float(val) == 1.0: is_pure = True
                except: pass
                
                if is_pure:
                    pure_points.append((row[lat_col], row[lon_col], idx))
        
        print(f"Found {len(pure_points)} pure points in SHP.")
        if len(pure_points) > 0:
            print("Sample pure points (Lat, Lon):")
            for p in pure_points[:10]:
                print(f"  ID {p[2]}: {p[0]}, {p[1]}")
                
    except Exception as e:
        print(f"Error loading SHP: {e}")
        return

    # 2. Check Downloads
    print("\nChecking Downloads...")
    if not os.path.exists(downloads_dir):
        print("No downloads dir")
        return
        
    folders = [f for f in os.listdir(downloads_dir) if os.path.isdir(os.path.join(downloads_dir, f))]
    print(f"Found {len(folders)} folders.")
    
    folder_coords = []
    for f in folders:
        parts = f.split('_')
        if parts[0] == 'L' and len(parts) > 3:
            try:
                lon = float(parts[1])
                lat = float(parts[2])
                folder_coords.append((lat, lon, f))
            except: pass
            
    print(f"Parsed {len(folder_coords)} folder coordinates.")
    
    # 3. Match
    print("\nAttempting Match...")
    match_count = 0
    pure_match_count = 0
    
    # We want to see closest matches for pure points
    for p_lat, p_lon, p_idx in pure_points:
        min_dist = 1e9
        best_match = None
        
        for f_lat, f_lon, f_name in folder_coords:
            dist = (p_lat - f_lat)**2 + (p_lon - f_lon)**2
            if dist < min_dist:
                min_dist = dist
                best_match = f_name
        
        print(f"Pure Point {p_idx} ({p_lat}, {p_lon}) closest match:")
        if best_match:
            print(f"  -> {best_match} (Dist sq: {min_dist:.10f})")
            if min_dist < (0.0001**2):
                print("  -> MATCHED!")
                pure_match_count += 1
            else:
                print("  -> NO MATCH (Too far)")
        else:
            print("  -> No folders found")

if __name__ == "__main__":
    diagnose()
