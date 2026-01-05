import geopandas as gpd
import pandas as pd

path = "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/training_punten_GEE_2.shp"
try:
    gdf = gpd.read_file(path)
    print(f"Total Rows: {len(gdf)}")
    print("Columns:", gdf.columns.tolist())
    
    # Print first 5 rows of all columns to see what we are dealing with
    print("\nFirst 5 rows:")
    print(gdf.head(5))
    
    # Specifically check the 'no soil' column if it exists
    cols = [c for c in gdf.columns if 'soil' in c.lower()]
    if cols:
        print(f"\nValue Counts for '{cols[0]}':")
        print(gdf[cols[0]].value_counts(dropna=False))
        
    # Check for anything looking like an ID
    id_candidates = [c for c in gdf.columns if 'id' in c.lower() or 'sample' in c.lower()]
    if id_candidates:
        print(f"\nID Candidates: {id_candidates}")
        for c in id_candidates:
            print(f"Head of {c}:")
            print(gdf[c].head(5))

except Exception as e:
    print(f"Error: {e}")
