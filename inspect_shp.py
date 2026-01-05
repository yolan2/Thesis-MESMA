import geopandas as gpd
import pandas as pd

SHP_PATH = "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/training_punten_GEE_2.shp"

try:
    gdf = gpd.read_file(SHP_PATH)
    print(f"Shapefile loaded. Rows: {len(gdf)}")
    print("Columns:", gdf.columns.tolist())
    
    # Find 'no_soil' column
    possible_names = ['no_soil', 'no soil', 'no.soil', 'no-soil', 'noSoil']
    col = next((c for c in gdf.columns if c.lower() in [p.lower() for p in possible_names]), None)
    
    if col:
        print(f"Found column: {col}")
        # Check indices for populus samples (e.g. 10, 11, 12 => indices 9, 10, 11)
        indices_to_check = [9, 10, 11, 12, 13, 14, 125, 126, 127, 128, 129]
        
        print(f"Checking specific indices (Sample ID = Index + 1):")
        for idx in indices_to_check:
            if idx < len(gdf):
                val = gdf.loc[idx, col]
                print(f"Index {idx} (Sample_{idx+1}): {col} = {val}")
            else:
                print(f"Index {idx} out of bounds")
    else:
        print("Column 'no_soil' not found.")
        
except Exception as e:
    print(f"Error: {e}")
