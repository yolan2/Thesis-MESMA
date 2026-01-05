import os
import geopandas as gpd
import pandas as pd

SHP_PATH = "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/training_punten_GEE_2.shp"
DOWNLOADS_DIR = "downloads"

def debug_data_loading():
    print("--- DEBUGGING DATA FILTERING ---")
    
    # 1. Load Shapefile and IDs
    if not os.path.exists(SHP_PATH):
        print(f"Shapefile not found: {SHP_PATH}")
        return

    gdf = gpd.read_file(SHP_PATH)
    print(f"Shapefile loaded. Rows: {len(gdf)}")
    
    # Identify ID column
    id_col = next((c for c in gdf.columns if c.lower() == 'id'), None)
    if id_col:
        print(f"Using ID column: {id_col}")
        shp_ids = set(gdf[id_col].astype(int))
    else:
        print("No ID column found. Using 1-based Index.")
        shp_ids = set(range(1, len(gdf) + 1))
        
    # Identify no_soil column
    possible_names = ['no_soil', 'no soil', 'no.soil', 'no-soil', 'noSoil']
    col = next((c for c in gdf.columns if c.lower() in [p.lower() for p in possible_names]), None)
    if not col:
        col = next((c for c in gdf.columns if 'no' in c.lower() and 'soil' in c.lower()), None)
    print(f"Using no_soil column: {col}")

    # Create mapping ID -> no_soil
    id_to_nosoil = {}
    if id_col:
        for _, row in gdf.iterrows():
            try:
                id_to_nosoil[int(row[id_col])] = row[col]
            except: pass
    else:
        for idx, row in gdf.iterrows():
            id_to_nosoil[idx + 1] = row[col]

    # 2. Iterate Folders
    if not os.path.exists(DOWNLOADS_DIR):
        print("Downloads dir not found.")
        return

    folders = [f for f in os.listdir(DOWNLOADS_DIR) if os.path.isdir(os.path.join(DOWNLOADS_DIR, f))]
    print(f"Found {len(folders)} folders.")
    
    counts = {
        "Total": 0,
        "Bad_Name_Format": 0,
        "Bad_Year": 0,
        "Type_None": 0,
        "Type_Excluded": 0,
        "ID_Not_In_SHP": 0,
        "No_Tiff": 0,
        "Valid_Pure": 0,
        "Valid_Impure": 0
    }
    
    sample_details = []

    for folder in folders:
        counts["Total"] += 1
        parts = folder.split('_')
        
        # 1. Parse ID
        try:
            sample_idx = int(parts[1])
        except:
            counts["Bad_Name_Format"] += 1
            continue
            
        # 2. Parse Type
        if 'SD' in parts:
            sd_index = parts.index('SD')
            typ = "_".join(parts[2:sd_index])
        elif len(parts) >= 3:
            typ = parts[2]
        else:
            counts["Bad_Name_Format"] += 1
            continue
            
        if typ == 'water': typ = 'barren'
        if typ.lower() == 'none':
            counts["Type_None"] += 1
            continue
        if typ in ['populus_herbs', 'tamarix_herbs']:
            counts["Type_Excluded"] += 1
            continue

        # 3. Parse Year
        year = None
        for p in parts:
            if len(p) == 4 and p.isdigit() and int(p) > 2000:
                year = int(p)
                break
        
        if not (year and 2023 <= year <= 2025):
            counts["Bad_Year"] += 1
            # print(f"Skipping {folder}: Bad Year {year}")
            continue

        # 4. Check SHP ID
        if sample_idx not in shp_ids:
            counts["ID_Not_In_SHP"] += 1
            print(f"Skipping {folder}: ID {sample_idx} not in SHP")
            continue
            
        # 5. Check Tiff
        sr_tif = None
        for root, dirs, files in os.walk(os.path.join(DOWNLOADS_DIR, folder)):
            for file in files:
                if 'AnalyticMS_SR' in file and file.endswith('.tif'):
                    sr_tif = os.path.join(root, file)
                    break
            if sr_tif: break
            
        if not sr_tif:
            counts["No_Tiff"] += 1
            # print(f"Skipping {folder}: No Tiff")
            continue
            
        # 6. Check Purity
        val = id_to_nosoil.get(sample_idx)
        is_pure = False
        if typ == 'barren':
            is_pure = True
        elif val == 1:
            is_pure = True
            
        if is_pure:
            counts["Valid_Pure"] += 1
            sample_details.append(f"{folder} (Pure, val={val})")
        else:
            counts["Valid_Impure"] += 1
            # sample_details.append(f"{folder} (Impure, val={val})")

    print("\n--- SUMMARY ---")
    for k, v in counts.items():
        print(f"{k}: {v}")
        
    print("\n--- PURE SAMPLES FOUND ---")
    # Group by type
    from collections import Counter
    types = [s.split('_')[2] for s in sample_details] # Rough parsing
    print(Counter(types))

if __name__ == "__main__":
    debug_data_loading()
