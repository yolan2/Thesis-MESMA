import os
import json
import time
import math
import requests
import geopandas as gpd
from shapely.geometry import Polygon
from pyproj import CRS, Transformer
from datetime import datetime, timedelta
from dateutil import parser
import zipfile
import rasterio

# --- 1. CREDENTIALS & PATHS ---
API_KEY = "PLAK5d5bab5be04c4cd38e620b159b991133"
GEE_PROJECT = "enhanced-oxygen-403615"
GEE_COLLECTION = "users/yolanrfclassifier/RF/planets"
INPUT_SHAPEFILE_PATH = "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/training_punten_GEE_2.shp"

# --- 2. SEARCH PARAMETERS ---
TARGET_DATES = ["2023-10-05T12:00:00Z", "2025-10-05T12:00:00Z"]
SEARCH_WINDOW_DAYS = 14     
CLOUD_COVER_MAX = 0.05      
MAX_CONCURRENT_ORDERS = 50

# --- 3. ENDPOINTS ---
SEARCH_URL = "https://api.planet.com/data/v1/quick-search"
ORDERS_URL = "https://api.planet.com/compute/ops/orders/v2"

def get_landsat_pixel_geom(lat, lon):
    """
    Calculates the 30m x 30m Landsat pixel footprint aligned to the UTM grid.
    """
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
    print(f"Columns in shapefile: {list(gdf.columns)}")
    if 'vegetation' in gdf.columns:
        unique_veg = gdf['vegetation'].unique()
        # Filter out None/NaN and sort, or convert to string to be safe
        safe_unique_veg = sorted([str(v) for v in unique_veg if v is not None])
        print(f"Vegetation types in shapefile: {safe_unique_veg}")
    else:
        print("No 'vegetation' column found in shapefile.")
    # Process all geometries in the geojson, regardless of vegetation type or no.soil
    # Removed filtering to include all points
    # Assume the geojson is in UTM zone 44N (based on coordinates around lon 82-88, lat 39-41)
    if gdf.crs is None:
        gdf.crs = 'EPSG:32644'
    # Reproject to WGS84
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
    
    print(f"Total Unique Pixels: {len(unique_geoms)}")
    return unique_geoms, unique_types

def search_planet_imagery(geom, target_date_str):
    target = parser.parse(target_date_str)
    date_config = {
        "gte": (target - timedelta(days=SEARCH_WINDOW_DAYS)).strftime("%Y-%m-%dT%H:%M:%S.%fZ"),
        "lte": (target + timedelta(days=SEARCH_WINDOW_DAYS)).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    }
    
    filter_params = {
        "type": "AndFilter",
        "config": [
            {"type": "GeometryFilter", "field_name": "geometry", "config": json.loads(json.dumps(geom.__geo_interface__))},
            {"type": "DateRangeFilter", "field_name": "acquired", "config": date_config},
            {"type": "RangeFilter", "field_name": "cloud_cover", "config": {"lte": CLOUD_COVER_MAX}},
            {"type": "StringInFilter", "field_name": "instrument", "config": ["PSB.SD"]}
        ]
    }

    res = requests.post(SEARCH_URL, auth=(API_KEY, ""), json={"item_types": ["PSScene"], "filter": filter_params})
    return res.json().get("features", []) if res.status_code == 200 else []

def get_closest_match(features, target_date_str):
    if not features: return []
    target = parser.parse(target_date_str)
    
    for f in features:
        f['time_diff'] = abs((parser.parse(f['properties']['acquired']) - target).total_seconds())
    
    sorted_features = sorted(features, key=lambda x: x['time_diff'])
    best_match = sorted_features[0]
    best_date_str = best_match['properties']['acquired'][0:10]
    days_diff = best_match['time_diff'] / 86400
    
    print(f"    -> Best Match: {best_date_str} (Off by {days_diff:.1f} days)")
    return [f['id'] for f in sorted_features if f['properties']['acquired'].startswith(best_date_str)]

def place_order_local(item_ids, geometry, name):
    order_payload = {
        "name": name,
        "products": [
            {
                "item_ids": item_ids,
                "item_type": "PSScene",
                "product_bundle": "analytic_8b_sr_udm2"
            }
        ],
        "delivery": {
            "single_archive": True,
            "archive_type": "zip"
        },
        "tools": [
            {"clip": {"aoi": json.loads(json.dumps(geometry.__geo_interface__))}}
        ]
    }

    headers = {'content-type': 'application/json'}
    
    # Retry loop for 429 Too Many Requests
    max_retries = 3
    for attempt in range(max_retries):
        res = requests.post(ORDERS_URL, auth=(API_KEY, ""), json=order_payload, headers=headers)
        
        if res.status_code == 202:
            order_id = res.json()['id']
            print(f"    [OK] Order Placed! ID: {order_id}")
            return order_id
        elif res.status_code == 429:
            print(f"    [WARN] Rate limit hit (429). Retrying in 5 seconds... (Attempt {attempt+1}/{max_retries})")
            time.sleep(5)
        else:
            print(f"    [ERROR] Order Failed: {res.text}")
            return None
    return None

def check_order_status(order_id):
    """
    Checks the status of a single order. Returns the JSON response if successful, or None on error.
    Does NOT block.
    """
    status_url = f"{ORDERS_URL}/{order_id}"
    res = requests.get(status_url, auth=(API_KEY, ""))
    if res.status_code != 200:
        print(f"    [ERROR] Failed to get order status for {order_id}: {res.text}")
        return None
    return res.json()

def download_order(order_details, download_dir="downloads"):
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)
    
    # Check if results exist
    if '_links' not in order_details or 'results' not in order_details['_links'] or not order_details['_links']['results']:
         print(f"    [ERROR] No download links found for order {order_details['id']}")
         return None

    download_url = order_details['_links']['results'][0]['location']
    filename = f"{order_details['name']}.zip"
    filepath = os.path.join(download_dir, filename)
    
    print(f"    Downloading to {filepath}...")
    res = requests.get(download_url, auth=(API_KEY, ""), stream=True)
    if res.status_code == 200:
        with open(filepath, 'wb') as f:
            for chunk in res.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"    [OK] Downloaded {filename}")
        
        # Unzip the file
        extract_dir = os.path.join(download_dir, order_details['name'])
        try:
            with zipfile.ZipFile(filepath, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)
            print(f"    [OK] Unzipped to {extract_dir}")
            return extract_dir
        except zipfile.BadZipFile:
             print(f"    [ERROR] Downloaded file is not a valid zip: {filepath}")
             return None
    else:
        print(f"    [ERROR] Download failed: {res.text}")
        return None

def export_to_gee(extract_dir, name):
    print(f"    Attempting to export {name} to GEE...")
    # For now, print instructions since direct upload requires cloud storage
    all_files = []
    for root, dirs, files in os.walk(extract_dir):
        for file in files:
            all_files.append(os.path.join(root, file))
    tif_files = [f for f in all_files if f.endswith('.tif')]
    if tif_files:
        print(f"    Found {len(tif_files)} TIFF files: {tif_files}")
    else:
        print(f"    No TIFF files found in {extract_dir}. Files present: {all_files}")

def monitor_active_orders(active_orders):
    """
    Checks statuses of active orders, downloads completed ones, and returns the updated list.
    """
    print(f"    [Monitor] Checking {len(active_orders)} active orders...")
    still_active = []
    
    for order in active_orders:
        order_id = order['id']
        name = order['name']
        
        details = check_order_status(order_id)
        if not details:
            # API error, keep checking? Or drop? Let's keep checking for transient errors.
            still_active.append(order)
            continue
            
        status = details['state']
        # print(f"    Order {name} ({order_id}): {status}") # Too verbose?
        
        if status == 'success':
            print(f"    -> Order {name} completed! Downloading...")
            download_path = download_order(details)
            if download_path:
                export_to_gee(download_path, name)
            # Done with this order
        elif status in ['failed', 'cancelled', 'partial']:
            print(f"    -> Order {name} failed/cancelled.")
            # Done with this order (failed)
        else:
            # queued, running, initializing -> keep waiting
            still_active.append(order)
    
    return still_active

import concurrent.futures

def process_search_task(task_args):
    """
    Worker function for parallel search.
    """
    i, geom, typ, target_date, download_dir = task_args
    year = target_date[:4]
    
    # Generate Location ID
    p = geom.centroid
    location_id = f"L_{p.x:.6f}_{p.y:.6f}"
    
    prefix = f"{location_id}_{typ}_SD_Oct5_{year}_"
    
    # Check existence and VALIDITY (8 bands)
    if os.path.exists(download_dir):
        existing_folders = [f for f in os.listdir(download_dir) if f.startswith(prefix)]
        for folder in existing_folders:
            folder_path = os.path.join(download_dir, folder)
            # Find the tif
            sr_tif = None
            for root, dirs, files in os.walk(folder_path):
                for file in files:
                    if 'AnalyticMS_SR' in file and file.endswith('.tif'):
                        sr_tif = os.path.join(root, file)
                        break
                if sr_tif: break
            
            if sr_tif:
                try:
                    with rasterio.open(sr_tif) as src:
                        if src.count == 8:
                            return {'status': 'exists', 'msg': f"{location_id} ({typ}) {year} already exists with 8 bands: {folder}"}
                        # else:
                        #     print(f"  Found {folder} but it has {src.count} bands. Re-downloading...")
                except Exception:
                    pass
            
    # Search
    try:
        # Rate limiting: simple sleep to avoid bursting too hard even with limited workers
        time.sleep(0.1) 
        features = search_planet_imagery(geom, target_date)
        if features:
            best_item_ids = get_closest_match(features, target_date)
            return {
                'status': 'found',
                'item_ids': best_item_ids,
                'geom': geom,
                'typ': typ,
                'year': year,
                'index': i,
                'location_id': location_id
            }
        else:
            return {'status': 'not_found', 'msg': f"No clear SuperDove imagery found for {location_id} ({typ}) in {year}."}
    except Exception as e:
        return {'status': 'error', 'msg': f"Error searching {location_id}: {str(e)}"}

def main():
    try:
        # Load all geometries and types
        geometries, types = prepare_pixel_aois(INPUT_SHAPEFILE_PATH)
    except Exception as e:
        print(f"Error: {e}")
        return

    print(f"\n--- STARTING DOWNLOAD RUN ---")
    
    if not geometries:
        print("Error: No geometries found in input file.")
        return

    print(f"Testing with {len(geometries)} sample(s).")

    active_orders = [] # List of dicts: {'id': order_id, 'name': order_name}
    download_dir = "downloads"
    
    # Prepare tasks
    search_tasks = []
    for target_date in TARGET_DATES:
        for i, (geom, typ) in enumerate(zip(geometries, types)):
            search_tasks.append((i, geom, typ, target_date, download_dir))

    print(f"\n--- PHASE 1: PARALLEL SEARCHING ({len(search_tasks)} tasks) ---")
    
    results_to_order = []
    
    # Use ThreadPoolExecutor to parallelize search (max_workers=5 to be safe with API limits)
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        # Submit all tasks
        future_to_task = {executor.submit(process_search_task, t): t for t in search_tasks}
        
        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_task):
            res = future.result()
            if res['status'] == 'found':
                print(f"    [Search Found] {res['location_id']} ({res['typ']}) {res['year']}")
                results_to_order.append(res)
            elif res['status'] == 'exists':
                print(f"    [Search Skip] {res['msg']}")
            elif res['status'] == 'not_found':
                # print(f"    [Search info] {res['msg']}") # Optional: reduce verbosity
                pass
            elif res['status'] == 'error':
                print(f"    [Search Error] {res['msg']}")

    print(f"\n--- PHASE 2: PLACING ORDERS ({len(results_to_order)} orders to place) ---")
    
    # Sort results by index/year just for consistent ordering if needed, though not strictly necessary
    results_to_order.sort(key=lambda x: (x['year'], x['index']))

    for res in results_to_order:
        # 1. Manage Concurrency
        while len(active_orders) >= MAX_CONCURRENT_ORDERS:
            print(f"    [Limit Reached] Waiting for some of the {len(active_orders)} active orders to complete...")
            active_orders = monitor_active_orders(active_orders)
            if len(active_orders) >= MAX_CONCURRENT_ORDERS:
                time.sleep(30) # Wait before re-checking
        
        # 2. Place Order
        i = res['index']
        typ = res['typ']
        year = res['year']
        geom = res['geom']
        item_ids = res['item_ids']
        location_id = res['location_id']
        
        order_name = f"{location_id}_{typ}_SD_Oct5_{year}_{int(time.time())}" 
        print(f"Placing order for {location_id} ({typ}) {year}...")
        order_id = place_order_local(item_ids, geom, order_name)
        
        if order_id:
            active_orders.append({'id': order_id, 'name': order_name})
            time.sleep(0.5) # Rate limit protection for Orders API
        else:
            print("    Order placement failed.")

    print(f"\n--- PHASE 3: FINISHING REMAINING {len(active_orders)} ORDERS ---")
    
    # --- PHASE 3: DRAIN REMAINING ORDERS ---
    while active_orders:
        active_orders = monitor_active_orders(active_orders)
        if active_orders:
            print(f"    Waiting 30 seconds before next check ({len(active_orders)} remaining)...")
            time.sleep(30)

    print(f"\n--- DOWNLOAD COMPLETE ---")

if __name__ == "__main__":
    main()
