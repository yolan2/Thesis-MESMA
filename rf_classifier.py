import os
import rasterio
from rasterio.warp import transform
import numpy as np
import pandas as pd
import geopandas as gpd
import re
from shapely.geometry import Point
from sklearn.ensemble import RandomForestClassifier, IsolationForest
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report
from collections import Counter

def parse_coordinate(val):
    """Parses coordinate values, handling '40,123' format and strings."""
    try:
        if pd.isna(val): return None
        if isinstance(val, (int, float)): return float(val)
        if isinstance(val, str):
            val = val.replace(',', '.')
            return float(val)
    except:
        return None
    return None

def collect_training_data(downloads_dir="downloads", spatial_reference_list=None):
    if spatial_reference_list is None: spatial_reference_list = []

    training_data = []
    labels = []
    groups = [] 
    group_purity = {} # Map group_id -> bool (True if pure) 
    
    sample_folders = [f for f in os.listdir(downloads_dir) if os.path.isdir(os.path.join(downloads_dir, f))]
    
    matched_count = 0
    total_folders = 0
    
    print(f"Scanning {len(sample_folders)} folders for training data...")
    
    # Pre-calculate pure point coordinates for debugging
    pure_coords = [p for p in spatial_reference_list if p[2] == 1]
    print(f"DEBUG: Tracking {len(pure_coords)} pure points from SHP for matching issues.")

    for folder in sample_folders:
        parts = folder.split('_')
        
        # Parse Coordinates from Folder Name (L_LON_LAT_TYPE...)
        folder_lat = None
        folder_lon = None
        typ = None
        
        is_l_format = (parts[0] == 'L' and len(parts) > 3)
        
        if is_l_format:
            try:
                # Format: L_{lon}_{lat}_{type}_...
                folder_lon = float(parts[1])
                folder_lat = float(parts[2])
                
                if 'SD' in parts:
                    sd_index = parts.index('SD')
                    typ = "_".join(parts[3:sd_index])
            except ValueError:
                continue
        else:
            # Fallback for Sample_ID formats, try to parse type at least
             if 'SD' in parts:
                sd_index = parts.index('SD')
                typ = "_".join(parts[2:sd_index])
             elif len(parts) >= 3:
                typ = parts[2]

        if not typ: continue
        
        # Standardize type
        if typ == 'water': typ = 'barren'
        if typ.lower() == 'none': continue

        allowed_classes = ['populus', 'barren', 'tamarix', 'herbs']
        
        # Check spatial match BEFORE type filtering for debug purposes
        is_pure = False
        matched_purity = 0
        match_found = False
        
        if folder_lat is not None and folder_lon is not None and spatial_reference_list:
            # Find closest point
            min_dist = 1e9
            
            for (ref_lat, ref_lon, ref_purity) in spatial_reference_list:
                dist = (folder_lat - ref_lat)**2 + (folder_lon - ref_lon)**2
                if dist < min_dist:
                    min_dist = dist
                    matched_purity = ref_purity
            
            # Threshold: 0.001 degrees is roughly 111 meters (relaxed from 0.00025)
            if min_dist < (0.001 ** 2):
                match_found = True
                is_pure = (matched_purity == 1)
        
        # Now apply type filtering
        if typ not in allowed_classes:
            # if match_found and is_pure:
            #     print(f"WARNING: Folder {folder} is spatially close to a PURE point but skipped due to type '{typ}'.")
            continue
            
        total_folders += 1
        
        if match_found:
             matched_count += 1
             if is_pure:
                 # print(f"  Match: {folder} -> Pure: {is_pure}")
                 pass
            
        group_id = folder
        group_purity[group_id] = is_pure
            
        group_id = folder
        group_purity[group_id] = is_pure

        # Find the analytic tif file recursively
        sr_tif = None
        for root, dirs, files in os.walk(os.path.join(downloads_dir, folder)):
            for file in files:
                if 'AnalyticMS_SR' in file and file.endswith('.tif'):
                    sr_tif = os.path.join(root, file)
                    break
            if sr_tif:
                break
        
        if not sr_tif:
            continue
            
        with rasterio.open(sr_tif) as src:
            if src.count != 8:
                continue
            
            bands = src.read()  # shape (8, h, w)
            mask = bands[0] != 0
            
            if not np.any(mask):
                continue
                
            valid_pixels = bands[:, mask].T 
            
            if valid_pixels.shape[0] > 0:
                training_data.extend(valid_pixels.tolist())
                labels.extend([typ] * valid_pixels.shape[0])
                groups.extend([group_id] * valid_pixels.shape[0])

    print(f"Matched {matched_count} folders to shapefile points out of {total_folders} candidates.")
    return np.array(training_data), np.array(labels), np.array(groups), group_purity

def remove_outliers(X, y, groups, contamination=0.05, random_state=42):
    """
    Removes outliers from the training data using Isolation Forest per class.
    """
    print("\n--- REMOVING OUTLIERS (Isolation Forest per Class) ---")
    print(f"Target Contamination: {contamination*100}%")
    
    unique_classes = np.unique(y)
    keep_mask = np.zeros(len(y), dtype=bool) 
    total_removed = 0
    
    for cls in unique_classes:
        cls_indices = np.where(y == cls)[0]
        X_cls = X[cls_indices]
        
        if len(X_cls) < 10:
            keep_mask[cls_indices] = True
            continue
            
        iso = IsolationForest(contamination=contamination, random_state=random_state, n_jobs=-1)
        preds = iso.fit_predict(X_cls)
        inliers_cls = (preds == 1)
        
        keep_indices = cls_indices[inliers_cls]
        keep_mask[keep_indices] = True
        
        n_removed = len(cls_indices) - len(keep_indices)
        total_removed += n_removed

    X_clean = X[keep_mask]
    y_clean = y[keep_mask]
    groups_clean = groups[keep_mask]
    
    print(f"Total Outliers Removed: {total_removed} / {len(y)} ({(total_removed/len(y))*100:.1f}%)")
    print(f"Cleaned Data Size: {len(X_clean)} pixels")
    return X_clean, y_clean, groups_clean

def purity_based_split(X, y, groups, group_purity):
    """
    Splits data into Training (Pure) and Testing (Impure).
    """
    unique_groups = np.unique(groups)
    train_groups = []
    test_groups = []
    
    for grp in unique_groups:
        if group_purity.get(grp, False):
            train_groups.append(grp)
        else:
            test_groups.append(grp)
            
    print(f"Split Strategy: {len(train_groups)} Pure Groups (Train) | {len(test_groups)} Impure Groups (Test)")
    
    train_mask = np.isin(groups, train_groups)
    test_mask = np.isin(groups, test_groups)
    return np.where(train_mask)[0], np.where(test_mask)[0]

def evaluate_rf_with_holdout(X, y, groups, group_purity):
    print("\n--- EVALUATING MODEL (Train on PURE, Test on IMPURE) ---")
    
    train_idx, test_idx = purity_based_split(X, y, groups, group_purity)
    
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    
    print(f"Training Pixels (Pure): {len(X_train)} | Test Pixels (Impure): {len(X_test)}")
    print(f"Training Class Distribution: {dict(Counter(y_train))}")
    
    if len(X_train) == 0:
        print("No pure training data. Using all data for training.")
        clf_final = RandomForestClassifier(n_estimators=100, max_depth=10, min_samples_split=10, min_samples_leaf=10, max_features='sqrt', class_weight='balanced', random_state=42, n_jobs=-1)
        clf_final.fit(X, y)
        return clf_final

    clf = RandomForestClassifier(n_estimators=100, max_depth=10, min_samples_split=10, min_samples_leaf=10, max_features='sqrt', class_weight='balanced', random_state=42, n_jobs=-1)
    clf.fit(X_train, y_train)
    
    if len(X_test) > 0:
        y_pred_test = clf.predict(X_test)
        std_acc = accuracy_score(y_test, y_pred_test)
        print(f"\n[Test] Standard Accuracy (Impure): {std_acc:.4f}")
        print("\nClassification Report (Standard Test):")
        print(classification_report(y_test, y_pred_test, labels=clf.classes_))
        print("\nConfusion Matrix (Standard Test):")
        print(confusion_matrix(y_test, y_pred_test, labels=clf.classes_))
        
        if 'barren' in clf.classes_:
            barren_idx = np.where(clf.classes_ == 'barren')[0][0]
            probs = clf.predict_proba(X_test)
            probs[:, barren_idx] = 0
            new_pred_indices = np.argmax(probs, axis=1)
            y_pred_veg = clf.classes_[new_pred_indices]
            veg_acc = accuracy_score(y_test, y_pred_veg)
            print(f"[Test] Vegetation-Only Accuracy (Impure): {veg_acc:.4f}")
            print("\nConfusion Matrix (Vegetation-Only Logic):")
            print(confusion_matrix(y_test, y_pred_veg, labels=clf.classes_))
            
    y_pred_train = clf.predict(X_train)
    train_acc = accuracy_score(y_train, y_pred_train)
    print(f"\n[Resubstitution] Training Pixel Accuracy (Pure): {train_acc:.4f}")
    print("\nConfusion Matrix (Training):")
    print(confusion_matrix(y_train, y_pred_train, labels=clf.classes_))

    clf_final = RandomForestClassifier(n_estimators=100, max_depth=10, min_samples_split=10, min_samples_leaf=10, max_features='sqrt', class_weight='balanced', random_state=42, n_jobs=-1)
    clf_final.fit(X_train, y_train)
    return clf_final

def predict_for_samples(downloads_dir, clf):
    print("\n--- PREDICTING FRACTIONAL COMPOSITION PER SAMPLE (Pixel-based) ---")
    
    sample_folders = [f for f in os.listdir(downloads_dir) if os.path.isdir(os.path.join(downloads_dir, f))]
    sample_folders.sort()
    
    predictions = []

    for folder in sample_folders:
        parts = folder.split('_')
        folder_lon = None
        folder_lat = None
        true_typ = None
        
        if parts[0] == 'L' and len(parts) > 3:
            try:
                folder_lon = float(parts[1])
                folder_lat = float(parts[2])
                if 'SD' in parts:
                    sd_index = parts.index('SD')
                    true_typ = "_".join(parts[3:sd_index])
            except:
                pass
        else:
             if 'SD' in parts:
                sd_index = parts.index('SD')
                true_typ = "_".join(parts[2:sd_index])
             elif len(parts) >= 3:
                true_typ = parts[2]

        year = None
        for p in parts:
            if len(p) == 4 and p.isdigit() and int(p) > 2000:
                year = int(p)
                break

        if true_typ == 'water': true_typ = 'barren'
        if not true_typ or true_typ.lower() == 'none': continue

        allowed_classes = ['populus', 'barren', 'tamarix', 'herbs']
        if true_typ not in allowed_classes: continue

        sr_tif = None
        for root, dirs, files in os.walk(os.path.join(downloads_dir, folder)):
            for file in files:
                if 'AnalyticMS_SR' in file and file.endswith('.tif'):
                    sr_tif = os.path.join(root, file)
                    break
            if sr_tif: break
        if not sr_tif: continue
            
        with rasterio.open(sr_tif) as src:
            if src.count != 8: continue
            
            bands = src.read()
            mask = bands[0] != 0
            if not np.any(mask): continue
                
            valid_pixels_data = bands[:, mask].T 
            preds = clf.predict(valid_pixels_data)
            counts = Counter(preds)
            total = len(preds)
            fractions = {k: v/total for k, v in counts.items()}
            
            loc_id = folder
            
            pred_row = {
                'location_id': loc_id,
                'year': year,
                'folder': folder,
                'true_type': true_typ
            }
            for cls in allowed_classes:
                pred_row[f'{cls}_frac'] = fractions.get(cls, 0.0)
            predictions.append(pred_row)
    
    if predictions:
        pred_df = pd.DataFrame(predictions)
        out_csv = "rf_predictions.csv"
        pred_df.to_csv(out_csv, index=False)
        print(f"\nSaved {len(predictions)} prediction records to {out_csv}")

def main():
    downloads_dir = "downloads"
    if not os.path.exists(downloads_dir): return

    shp_path_candidates = [
        "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/training_punten_GEE_3.shp",
        "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/training_punten_GEE_2.shp",
    ]
    shp_path = next((p for p in shp_path_candidates if os.path.exists(p)), shp_path_candidates[0])
    spatial_reference_list = []
    
    if os.path.exists(shp_path):
        print(f"Loading shapefile from {shp_path}...")
        try:
            gdf = gpd.read_file(shp_path)
            gdf.columns = [c.lower() for c in gdf.columns]
            print(f"Shapefile Columns: {gdf.columns.tolist()}")
            
            # Identify columns
            soil_col = next((c for c in gdf.columns if 'soil' in c or 'pure' in c), None)

            # Some GEE exports include lat/lon attributes even when geometry is NULL.
            lat_col = next((c for c in gdf.columns if c in ("lat", "latitude", "y", "lat_dd", "lat_deg", "lat_wgs84")), None)
            lon_col = next((c for c in gdf.columns if c in ("lon", "lng", "longitude", "x", "lon_dd", "lon_deg", "lon_wgs84")), None)
            if lat_col and lon_col:
                print(f"Using fallback coordinates from columns: '{lat_col}', '{lon_col}'")
            
            if soil_col:
                print(f"Using Purity Column: '{soil_col}'")
                print(f"Total rows in Shapefile: {len(gdf)}")
                print(f"Value Counts for '{soil_col}':\n{gdf[soil_col].value_counts(dropna=False)}")
                
                # Ensure we are in WGS84 to match the L_LON_LAT format
                if gdf.crs and gdf.crs.to_string() != "EPSG:4326":
                    print(f"Reprojecting from {gdf.crs.to_string()} to EPSG:4326...")
                    gdf = gdf.to_crs(epsg=4326)
                
                count = 0
                dropped_stats = {"missing_geom": 0, "parse_error": 0, "filled_from_latlon": 0, "missing_latlon": 0}
                
                for idx, row in gdf.iterrows():
                    try:
                        # Prefer geometry, but if it's missing use lat/lon attribute columns (GEE_3 export).
                        lat = None
                        lon = None

                        if row.geometry is None or row.geometry.is_empty:
                            if lat_col and lon_col:
                                lat = parse_coordinate(row.get(lat_col))
                                lon = parse_coordinate(row.get(lon_col))
                                if lat is not None and lon is not None:
                                    dropped_stats["filled_from_latlon"] += 1
                                else:
                                    dropped_stats["missing_latlon"] += 1
                                    dropped_stats["missing_geom"] += 1
                                    continue
                            else:
                                dropped_stats["missing_geom"] += 1
                                continue
                        else:
                            # Use centroid for robustness if not Point (though likely Point)
                            centroid = row.geometry.centroid
                            lat = centroid.y
                            lon = centroid.x
                        
                        raw_soil = row[soil_col]
                        # Handle 1.0 vs NaN
                        is_pure = 0
                        if pd.notna(raw_soil):
                            try:
                                # Check for various forms of "1" (1, 1.0, '1', 100, etc)
                                val_float = float(raw_soil)
                                if val_float == 1.0 or val_float == 100.0: is_pure = 1
                            except:
                                pass
                        
                        spatial_reference_list.append((lat, lon, is_pure))
                        count += 1
                    except Exception as e:
                        dropped_stats["parse_error"] += 1
                        pass
                
                print(f"Loaded {len(spatial_reference_list)} reference points.")
                print(f"Pure points (val=1): {sum(p[2] for p in spatial_reference_list)}")
                
                if len(spatial_reference_list) < len(gdf):
                    print(f"Dropped {len(gdf) - len(spatial_reference_list)} rows. Reasons: {dropped_stats}")
            else:
                print("Could not find purity column (no_soil) in shapefile.")
        except Exception as e:
            print(f"Error loading shapefile: {e}")
    
    print("\n--- COLLECTING DATA ---")
    # Threshold: 0.001 degrees is roughly 111 meters
    X, y, groups, group_purity = collect_training_data(downloads_dir, spatial_reference_list)

    # --- DIAGNOSTIC: Check why pure points didn't match ---
    # pure_refs = [p for p in spatial_reference_list if p[2] == 1]
    # if pure_refs:
    #     print(f"\n--- DIAGNOSING {len(pure_refs)} PURE POINTS MATCHING ---")
    #     sample_folders = [f for f in os.listdir(downloads_dir) if os.path.isdir(os.path.join(downloads_dir, f))]
    #     
    #     for i, (plat, plon, _) in enumerate(pure_refs):
    #         min_dist = 1e9
    #         best_folder = None
    #         
    #         for folder in sample_folders:
    #             parts = folder.split('_')
    #             if parts[0] == 'L' and len(parts) > 3:
    #                 try:
    #                     flon = float(parts[1])
    #                     flat = float(parts[2])
    #                     dist = (plat - flat)**2 + (plon - flon)**2
    #                     if dist < min_dist:
    #                         min_dist = dist
    #                         best_folder = folder
    #                 except: pass
    #         
    #         dist_deg = np.sqrt(min_dist)
    #         print(f"Pure Point {i+1} ({plat:.5f}, {plon:.5f}) closest: {best_folder}")
    #         print(f"   -> Dist: {dist_deg:.6f} deg")
    #         if dist_deg > 0.00025:
    #             print("   -> FAIL: Too far (Threshold 0.00025)")
    #         else:
    #             print("   -> OK: Within threshold (Check type filtering?)")
    # # ------------------------------------------------------
    
    # Special handling for barren: split barren groups in half for pure/impure
    barren_groups = [grp for grp in np.unique(groups) if any(y[i] == 'barren' for i in np.where(groups == grp)[0])]
    if barren_groups:
        print(f"Barren groups found: {len(barren_groups)}. Splitting in half for pure/impure.")
        seed = int(os.environ.get("MESMA_SEED", "42"))
        np.random.seed(seed)
        import random as _py_random
        _py_random.seed(seed)
        np.random.shuffle(barren_groups)
        half = len(barren_groups) // 2
        pure_barren = barren_groups[:half]
        impure_barren = barren_groups[half:]
        for grp in pure_barren:
            group_purity[grp] = True
        for grp in impure_barren:
            group_purity[grp] = False
        print(f"Pure barren groups: {len(pure_barren)}, Impure barren groups: {len(impure_barren)}")

    if len(X) == 0:
        print("No training data found.")
        return
        
    X, y, groups = remove_outliers(X, y, groups, contamination=0.05)
    
    if len(np.unique(groups)) < 1: return
    
    clf = evaluate_rf_with_holdout(X, y, groups, group_purity)
    print("Final model trained.")

    predict_for_samples(downloads_dir, clf)
    print("\n--- PROCESS COMPLETE ---")

if __name__ == "__main__":
    main()
