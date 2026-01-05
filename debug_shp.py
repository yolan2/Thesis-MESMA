import geopandas as gpd
import pandas as pd

SHP_PATH = "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/training_punten_GEE_2.shp"

try:
    gdf = gpd.read_file(SHP_PATH)
    print("Columns:", gdf.columns.tolist())
    print(f"Shape: {gdf.shape}")
    
    # Check for ID columns
    id_cols = [c for c in gdf.columns if 'id' in c.lower()]
    print("Potential ID columns:", id_cols)
    
    # Check no_soil column
    possible_names = ['no_soil', 'no soil', 'no.soil', 'no-soil', 'noSoil']
    no_soil_col = next((c for c in gdf.columns if c.lower() in [p.lower() for p in possible_names]), None)
    
    if no_soil_col:
        print(f"Found 'no_soil' column: {no_soil_col}")
        print(gdf[no_soil_col].value_counts(dropna=False))
        
        # Check first 5 rows
        print("\nFirst 5 rows:")
        cols_to_show = id_cols + [no_soil_col]
        print(gdf[cols_to_show].head(10))
        
        # Check specific samples mentioned in previous logs (e.g. Sample 208 => Index 207?)
        # Logs said: "Sample_208... (Idx 207) ... value is 'nan'"
        print("\nChecking Index 207:")
        if 207 in gdf.index:
            print(gdf.loc[207, cols_to_show])
        else:
            print("Index 207 not found.")

    else:
        print("'no_soil' column not found.")

except Exception as e:
    print(f"Error: {e}")
