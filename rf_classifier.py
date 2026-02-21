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

# --- Align RF feature space with MESMA (RAW_BANDS + OPTIMAL_INDICES) ---
MESMA_RAW_BANDS = ["blue", "green", "red", "nir", "swir1", "swir2"]
MESMA_OPTIMAL_INDICES = ["PSRI", "NDMI", "EVI", "TCW", "NDTI", "BSI", "MSI", "SIPI"]


def _compute_mesma_indices(band_df: "pd.DataFrame") -> "pd.DataFrame":
    """Compute MESMA-style indices from available band columns in `band_df`.

    - `band_df` should contain one or more of the columns named in MESMA_RAW_BANDS.
    - Returns a DataFrame with any available RAW_BANDS followed by the MESMA_OPTIMAL_INDICES
      that can be computed from the available bands.
    - Uses the same formulas as `compute_indices_from_bands()` in `fit_veg_mixture_mesma.R`.
    """
    eps = 1e-9
    out = pd.DataFrame(index=band_df.index)

    # copy available raw bands (keep same order as MESMA_RAW_BANDS)
    for b in MESMA_RAW_BANDS:
        if b in band_df.columns:
            out[b] = pd.to_numeric(band_df[b], errors="coerce")

    # convenience local alias access
    b = out

    # PSRI: (red - blue) / (nir + eps)
    if set(["red", "blue", "nir"]).issubset(b.columns):
        out["PSRI"] = (b["red"] - b["blue"]) / (b["nir"].astype(float) + eps)

    # NDMI: (nir - swir1) / (nir + swir1 + eps)
    if set(["nir", "swir1"]).issubset(b.columns):
        out["NDMI"] = (b["nir"] - b["swir1"]) / (b["nir"] + b["swir1"] + eps)

    # EVI: 2.5 * ((nir - red) / (nir + 6*red - 7.5*blue + 1 + eps))
    if set(["nir", "red", "blue"]).issubset(b.columns):
        out["EVI"] = 2.5 * ((b["nir"] - b["red"]) / (b["nir"] + 6*b["red"] - 7.5*b["blue"] + 1 + eps))

    # TCW: (swir1 - swir2) / (swir1 + swir2 + eps)
    if set(["swir1", "swir2"]).issubset(b.columns):
        out["TCW"] = (b["swir1"] - b["swir2"]) / (b["swir1"] + b["swir2"] + eps)

    # NDTI: (swir1 - swir2) / (swir1 + swir2 + eps)  (same form as TCW)
    if set(["swir1", "swir2"]).issubset(b.columns):
        out["NDTI"] = (b["swir1"] - b["swir2"]) / (b["swir1"] + b["swir2"] + eps)

    # BSI: ((swir1 + red) - (nir + blue)) / ((swir1 + red) + (nir + blue) + eps)
    if set(["swir1", "red", "nir", "blue"]).issubset(b.columns):
        term1 = b["swir1"] + b["red"]
        term2 = b["nir"] + b["blue"]
        out["BSI"] = (term1 - term2) / (term1 + term2 + eps)

    # MSI: swir1 / (nir + eps)
    if set(["swir1", "nir"]).issubset(b.columns):
        out["MSI"] = b["swir1"] / (b["nir"] + eps)

    # SIPI: (nir - blue) / (nir - red + eps)
    if set(["nir", "blue", "red"]).issubset(b.columns):
        out["SIPI"] = (b["nir"] - b["blue"]) / (b["nir"] - b["red"] + eps)

    # Keep column order predictable: raw bands (in MESMA_RAW_BANDS order) then MESMA_OPTIMAL_INDICES
    cols = [c for c in MESMA_RAW_BANDS if c in out.columns]
    cols += [c for c in MESMA_OPTIMAL_INDICES if c in out.columns]
    return out.loc[:, cols]


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

    rows = []           # list of dicts (feature-name -> value)
    labels = []
    groups = []
    group_purity = {}

    sample_folders = [f for f in os.listdir(downloads_dir) if os.path.isdir(os.path.join(downloads_dir, f))]
    print(f"Scanning {len(sample_folders)} folders for training data...")

    # Pre-calc pure points for diagnostics
    pure_coords = [p for p in spatial_reference_list if p[2] == 1]
    print(f"DEBUG: Tracking {len(pure_coords)} pure points from SHP for matching issues.")

    for folder in sample_folders:
        parts = folder.split('_')
        # parse folder name (same logic as before)
        folder_lat = None; folder_lon = None; typ = None
        is_l_format = (parts[0] == 'L' and len(parts) > 3)
        if is_l_format:
            try:
                folder_lon = float(parts[1]); folder_lat = float(parts[2])
                if 'SD' in parts:
                    sd_index = parts.index('SD'); typ = "_".join(parts[3:sd_index])
            except ValueError:
                continue
        else:
            if 'SD' in parts:
                sd_index = parts.index('SD'); typ = "_".join(parts[2:sd_index])
            elif len(parts) >= 3:
                typ = parts[2]
        if not typ: continue
        if typ == 'water': typ = 'barren'
        if typ.lower() == 'none': continue

        allowed_classes = ['populus', 'barren', 'tamarix', 'herbs']
        # spatial match diagnostic
        is_pure = False; match_found = False; matched_purity = 0
        if folder_lat is not None and folder_lon is not None and spatial_reference_list:
            min_dist = 1e9
            for (ref_lat, ref_lon, ref_purity) in spatial_reference_list:
                dist = (folder_lat - ref_lat)**2 + (folder_lon - ref_lon)**2
                if dist < min_dist:
                    min_dist = dist; matched_purity = ref_purity
            if min_dist < (0.001 ** 2):
                match_found = True; is_pure = (matched_purity == 1)
        if typ not in allowed_classes: continue

        group_id = folder
        group_purity[group_id] = is_pure

        # find analytic tif
        sr_tif = None
        for root, dirs, files in os.walk(os.path.join(downloads_dir, folder)):
            for file in files:
                if 'AnalyticMS_SR' in file and file.endswith('.tif'):
                    sr_tif = os.path.join(root, file); break
            if sr_tif: break
        if not sr_tif: continue

        with rasterio.open(sr_tif) as src:
            # read bands and try to construct named-band DataFrame when possible
            bands = src.read()  # shape (nb, h, w)
            mask = bands[0] != 0
            if not np.any(mask): continue
            valid_pixels = bands[:, mask].T  # (n_pixels, n_bands)

            # attempt to detect band -> name mapping
            band_names = []
            descs = []
            try:
                descs = [d.lower() if d else '' for d in (src.descriptions or [])]
            except Exception:
                descs = []

            if descs and len(descs) == src.count:
                # try to map known names from descriptions
                for d in descs:
                    if 'blue' in d: band_names.append('blue')
                    elif 'green' in d: band_names.append('green')
                    elif 'red' in d and 'red edge' not in d: band_names.append('red')
                    elif 'nir' in d and 'swir' not in d: band_names.append('nir')
                    elif 'swir1' in d or 'swir 1' in d or 'swir' in d and 'swir2' not in d: band_names.append('swir1')
                    elif 'swir2' in d or 'swir 2' in d: band_names.append('swir2')
                    else:
                        band_names.append(None)
            else:
                # heuristic: if image has exactly 6 bands assume Landsat-like order
                if src.count == 6:
                    band_names = ['blue','green','red','nir','swir1','swir2']
                else:
                    # fallback: create generic band names (preserve ordering)
                    band_names = [f'band_{i+1}' for i in range(src.count)]

            # build DataFrame for the raw bands (columns named)
            df_bands = pd.DataFrame(valid_pixels, columns=band_names)

            # try to compute MESMA indices where possible; if not possible we'll keep raw bands
            fea_df = _compute_mesma_indices(df_bands)
            if fea_df.shape[1] == 0:
                # fallback to raw bands (use generic names)
                fea_df = df_bands.copy()

            # append rows
            for i, row in fea_df.iterrows():
                rows.append(row.to_dict()); labels.append(typ); groups.append(group_id)

    print(f"Collected features from {len(rows)} pixels across {len(set(groups))} groups (folders)")
    if len(rows) == 0:
        return pd.DataFrame(), np.array([]), np.array([]), group_purity

    X_df = pd.DataFrame(rows)
    # Ensure deterministic column order: MESMA raw bands first then OPTIMAL_INDICES then any remaining
    cols_order = [c for c in MESMA_RAW_BANDS if c in X_df.columns]
    cols_order += [c for c in MESMA_OPTIMAL_INDICES if c in X_df.columns and c not in cols_order]
    # add any other columns at end
    cols_order += [c for c in X_df.columns if c not in cols_order]
    X_df = X_df.loc[:, cols_order]

    # diagnostic: print final feature list so user can confirm consistency with MESMA
    print(f"[FEATURE SPACE] Training features: {list(X_df.columns)}")

    return X_df, np.array(labels), np.array(groups), group_purity

def remove_outliers(X, y, groups, contamination=0.05, random_state=42):
    """
    Removes outliers from the training data using Isolation Forest per class.

    Accepts X as a numpy array or pandas DataFrame. Returns X in the same *type*
    as the input (DataFrame preserved when given).
    """
    print("\n--- REMOVING OUTLIERS (Isolation Forest per Class) ---")
    print(f"Target Contamination: {contamination*100}%")

    is_df = isinstance(X, pd.DataFrame)
    n_rows = len(y)
    unique_classes = np.unique(y)
    keep_mask = np.zeros(n_rows, dtype=bool)
    total_removed = 0

    for cls in unique_classes:
        cls_indices = np.where(y == cls)[0]
        if is_df:
            X_cls = X.iloc[cls_indices]
        else:
            X_cls = X[cls_indices]

        if len(X_cls) < 10:
            keep_mask[cls_indices] = True
            continue

        iso = IsolationForest(contamination=contamination, random_state=random_state, n_jobs=-1)
        # sklearn accepts DataFrame or ndarray
        preds = iso.fit_predict(X_cls)
        inliers_cls = (preds == 1)

        keep_indices = cls_indices[inliers_cls]
        keep_mask[keep_indices] = True

        n_removed = len(cls_indices) - len(keep_indices)
        total_removed += n_removed

    if is_df:
        X_clean = X.iloc[keep_mask].reset_index(drop=True)
    else:
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

    # Accept X as DataFrame or ndarray
    is_df = isinstance(X, pd.DataFrame)
    if is_df:
        X_df = X.copy().reset_index(drop=True)
    else:
        X_df = pd.DataFrame(X)

    train_idx, test_idx = purity_based_split(X_df.values, y, groups, group_purity)

    X_train_df = X_df.iloc[train_idx]
    X_test_df = X_df.iloc[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    print(f"Training Pixels (Pure): {len(X_train_df)} | Test Pixels (Impure): {len(X_test_df)}")
    print(f"Training Class Distribution: {dict(Counter(y_train))}")

    if len(X_train_df) == 0:
        print("No pure training data. Using all data for training.")
        clf_final = RandomForestClassifier(n_estimators=100, max_depth=10, min_samples_split=10, min_samples_leaf=10, max_features='sqrt', class_weight='balanced', random_state=42, n_jobs=-1)
        # pass DataFrame so sklearn records feature names
        clf_final.fit(X_df, y)
        return clf_final

    clf = RandomForestClassifier(n_estimators=100, max_depth=10, min_samples_split=10, min_samples_leaf=10, max_features='sqrt', class_weight='balanced', random_state=42, n_jobs=-1)
    clf.fit(X_train_df, y_train)

    if len(X_test_df) > 0:
        y_pred_test = clf.predict(X_test_df)
        std_acc = accuracy_score(y_test, y_pred_test)
        print(f"\n[Test] Standard Accuracy (Impure): {std_acc:.4f}")
        print("\nClassification Report (Standard Test):")
        print(classification_report(y_test, y_pred_test, labels=clf.classes_))
        print("\nConfusion Matrix (Standard Test):")
        print(confusion_matrix(y_test, y_pred_test, labels=clf.classes_))

        if 'barren' in clf.classes_:
            barren_idx = np.where(clf.classes_ == 'barren')[0][0]
            probs = clf.predict_proba(X_test_df)
            probs[:, barren_idx] = 0
            new_pred_indices = np.argmax(probs, axis=1)
            y_pred_veg = clf.classes_[new_pred_indices]
            veg_acc = accuracy_score(y_test, y_pred_veg)
            print(f"[Test] Vegetation-Only Accuracy (Impure): {veg_acc:.4f}")
            print("\nConfusion Matrix (Vegetation-Only Logic):")
            print(confusion_matrix(y_test, y_pred_veg, labels=clf.classes_))

    y_pred_train = clf.predict(X_train_df)
    train_acc = accuracy_score(y_train, y_pred_train)
    print(f"\n[Resubstitution] Training Pixel Accuracy (Pure): {train_acc:.4f}")
    print("\nConfusion Matrix (Training):")
    print(confusion_matrix(y_train, y_pred_train, labels=clf.classes_))

    clf_final = RandomForestClassifier(n_estimators=100, max_depth=10, min_samples_split=10, min_samples_leaf=10, max_features='sqrt', class_weight='balanced', random_state=42, n_jobs=-1)
    clf_final.fit(X_train_df, y_train)
def predict_for_samples(downloads_dir, clf):
    print("\n--- PREDICTING FRACTIONAL COMPOSITION PER SAMPLE (Pixel-based) ---")

    # Attempt to retrieve training feature order from the trained classifier
    feature_cols = None
    if hasattr(clf, 'feature_names_in_'):
        try:
            feature_cols = list(clf.feature_names_in_)
        except Exception:
            feature_cols = None

    sample_folders = [f for f in os.listdir(downloads_dir) if os.path.isdir(os.path.join(downloads_dir, f))]
    sample_folders.sort()

    predictions = []

    for folder in sample_folders:
        parts = folder.split('_')
        folder_lon = None; folder_lat = None; true_typ = None
        if parts[0] == 'L' and len(parts) > 3:
            try:
                folder_lon = float(parts[1]); folder_lat = float(parts[2])
                if 'SD' in parts:
                    sd_index = parts.index('SD'); true_typ = "_".join(parts[3:sd_index])
            except:
                pass
        else:
            if 'SD' in parts:
                sd_index = parts.index('SD'); true_typ = "_".join(parts[2:sd_index])
            elif len(parts) >= 3:
                true_typ = parts[2]

        year = None
        for p in parts:
            if len(p) == 4 and p.isdigit() and int(p) > 2000:
                year = int(p); break

        if true_typ == 'water': true_typ = 'barren'
        if not true_typ or true_typ.lower() == 'none': continue

        allowed_classes = ['populus', 'barren', 'tamarix', 'herbs']
        if true_typ not in allowed_classes: continue

        sr_tif = None
        for root, dirs, files in os.walk(os.path.join(downloads_dir, folder)):
            for file in files:
                if 'AnalyticMS_SR' in file and file.endswith('.tif'):
                    sr_tif = os.path.join(root, file); break
            if sr_tif: break
        if not sr_tif: continue

        with rasterio.open(sr_tif) as src:
            bands = src.read()
            mask = bands[0] != 0
            if not np.any(mask): continue
            valid_pixels = bands[:, mask].T  # (n_pixels, n_bands)

            # attempt to name bands like during training
            descs = []
            try:
                descs = [d.lower() if d else '' for d in (src.descriptions or [])]
            except Exception:
                descs = []
            if descs and len(descs) == src.count:
                band_names = []
                for d in descs:
                    if 'blue' in d: band_names.append('blue')
                    elif 'green' in d: band_names.append('green')
                    elif 'red' in d and 'red edge' not in d: band_names.append('red')
                    elif 'nir' in d and 'swir' not in d: band_names.append('nir')
                    elif 'swir1' in d or 'swir 1' in d or ('swir' in d and 'swir2' not in d): band_names.append('swir1')
                    elif 'swir2' in d or 'swir 2' in d: band_names.append('swir2')
                    else: band_names.append(None)
            else:
                if src.count == 6:
                    band_names = ['blue','green','red','nir','swir1','swir2']
                else:
                    band_names = [f'band_{i+1}' for i in range(src.count)]

            df_bands = pd.DataFrame(valid_pixels, columns=band_names)
            fea_df = _compute_mesma_indices(df_bands)
            if fea_df.shape[1] == 0:
                fea_df = df_bands.copy()

            # Align with training feature columns if available
            if feature_cols is not None:
                # ensure all required cols present; fill missing with 0
                missing = [c for c in feature_cols if c not in fea_df.columns]
                if missing:
                    # safer to fill with zeros rather than drop
                    for m in missing:
                        fea_df[m] = 0.0
                fea_df = fea_df.reindex(columns=feature_cols)

            # predict (sklearn accepts DataFrame directly)
            try:
                preds = clf.predict(fea_df)
            except Exception:
                # fallback: convert to numpy array
                preds = clf.predict(fea_df.values)

            counts = Counter(preds)
            total = len(preds)
            fractions = {k: v/total for k, v in counts.items()}

            loc_id = folder
            pred_row = {'location_id': loc_id, 'year': year, 'folder': folder, 'true_type': true_typ}
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
        seed = int(os.environ.get("MESMA_SEED", "123"))
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
