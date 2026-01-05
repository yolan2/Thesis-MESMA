import os
import geopandas as gpd
from shapely.geometry import Polygon
from pyproj import CRS, Transformer
import math

INPUT_SHAPEFILE_PATH = "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/training_punten_GEE_2.shp"
DOWNLOAD_DIR = "downloads"

def get_landsat_pixel_geom(lat, lon):
    # 1. Determine UTM Zone
    zone_number = int((lon + 180) / 6) + 1
    is_northern = lat >= 0
    
    # 2. Define Projections
    wgs84_crs = CRS.from_epsg(4326)
    hemisphere = "north" if is_northern else "south"
    utm_crs_str = f"+proj=utm +zone={zone_number} +{hemisphere} +datum=WGS84 +units=m +no_defs"
    utm_crs = CRS.from_string(utm_crs_str)
    
    # Create Transformers
    to_utm = Transformer.from_crs(wgs84_crs, utm_crs, always_xy=True)
    to_wgs84 = Transformer.from_crs(utm_crs, wgs84_crs, always_xy=True)
    
    # 3. Project Point to UTM
    easting, northing = to_utm.transform(lon, lat)
    
    # 4. Snap to 30m Landsat Grid
    e_min = math.floor(easting / 30) * 30
    n_min = math.floor(northing / 30) * 30
    e_max = e_min + 30
    n_max = n_min + 30
    
    # 5. Create 30x30m Square
    utm_poly_coords = [
        (e_min, n_min), (e_max, n_min), 
        (e_max, n_max), (e_min, n_max), 
        (e_min, n_min)
    ]
    
    # 6. Reproject back to WGS84
    wgs84_coords = [to_wgs84.transform(x, y) for x, y in utm_poly_coords]
    return Polygon(wgs84_coords)

def prepare_pixel_aois(input_shapefile_path):
    print(f"--- GENERATING LANDSAT PIXEL ALIGNMENTS ---")
    if not os.path.exists(input_shapefile_path):
        raise FileNotFoundError(f"{input_shapefile_path} not found.")
        
    gdf = gpd.read_file(input_shapefile_path)
    if gdf.crs is None:
        gdf.crs = 'EPSG:32644'
    gdf = gdf.to_crs('EPSG:4326')
    
    pixel_data = []
    for idx, row in gdf.iterrows():
        if row.geometry is None or row.geometry.is_empty:
            continue
        p = row.geometry.centroid
        pixel_geom = get_landsat_pixel_geom(p.y, p.x)
        pixel_data.append({'geom': pixel_geom, 'type': row.get('vegetation', 'unknown')})
    
    pixel_gdf = gpd.GeoDataFrame(pixel_data)
    
    # Remove duplicates based on geometry
    pixel_gdf['wkt'] = pixel_gdf['geom'].apply(lambda g: g.wkt)
    unique_pixel_data = pixel_gdf.drop_duplicates(subset=['wkt'])
    
    unique_geoms = unique_pixel_data['geom'].tolist()
    unique_types = unique_pixel_data['type'].tolist()
    
    return unique_geoms, unique_types

def main():
    try:
        geometries, types = prepare_pixel_aois(INPUT_SHAPEFILE_PATH)
    except Exception as e:
        print(f"Error loading shapefile: {e}")
        return

    print(f"Found {len(geometries)} unique geometries.")
    
    # Create mapping: Index -> Lat/Lon ID
    mapping = {}
    for i, geom in enumerate(geometries):
        centroid = geom.centroid
        # Format: L_LON_LAT
        lon_str = f"{centroid.x:.6f}"
        lat_str = f"{centroid.y:.6f}"
        new_id = f"L_{lon_str}_{lat_str}"
        mapping[i + 1] = new_id  # Sample_1 is index 0
        
        # print(f"Sample {i+1} -> {new_id}")

    # Process downloads directory
    if not os.path.exists(DOWNLOAD_DIR):
        print(f"{DOWNLOAD_DIR} does not exist.")
        return

    files = os.listdir(DOWNLOAD_DIR)
    renamed_count = 0
    
    # Sort files to ensure we don't process renamed files or create conflicts easily
    # But renaming "Sample_1_..." won't conflict with "L_..."
    
    for filename in files:
        if not filename.startswith("Sample_"):
            continue
            
        # Parse Sample ID
        try:
            parts = filename.split('_')
            # Sample_1_barren...
            # parts[0] = "Sample"
            # parts[1] = "1"
            sample_idx = int(parts[1])
        except ValueError:
            print(f"Skipping malformed filename: {filename}")
            continue
            
        if sample_idx in mapping:
            new_id = mapping[sample_idx]
            # Reconstruct remainder
            # Sample_1_... -> L_..._...
            remainder = "_".join(parts[2:])
            new_filename = f"{new_id}_{remainder}"
            
            old_path = os.path.join(DOWNLOAD_DIR, filename)
            new_path = os.path.join(DOWNLOAD_DIR, new_filename)
            
            # Rename
            print(f"Renaming {filename} -> {new_filename}")
            try:
                os.rename(old_path, new_path)
                renamed_count += 1
            except OSError as e:
                print(f"Error renaming {filename}: {e}")
        else:
            print(f"Warning: Sample index {sample_idx} not found in shapefile mapping.")

    print(f"Renamed {renamed_count} files/folders.")

if __name__ == "__main__":
    main()
