#!/usr/bin/env python3
import os
import glob
import argparse
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from xgboost import XGBRegressor


def load_labels(csv_path: str):
    df = pd.read_csv(csv_path)
    return dict(zip(df['id'], df['hc50_ug_ml']))


def load_feature_set(features_root: str):
    rows = []
    for sample_dir in glob.glob(os.path.join(features_root, '*')):
        sample_id = os.path.basename(sample_dir)
        reps = sorted(glob.glob(os.path.join(sample_dir, 'rep*', 'features.csv')))
        if not reps:
            continue
        # Mittelung über letzte Fensterwerte (bereits Fenster-beschränkt)
        vals = []
        for rep_file in reps:
            df = pd.read_csv(rep_file)
            vals.append(df[['insertion_nm','contacts_6A']].mean().values)
        mean_vals = np.mean(np.stack(vals, axis=0), axis=0)
        rows.append({'id': sample_id, 'insertion_nm': mean_vals[0], 'contacts_6A': mean_vals[1]})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels_csv', default='data/raw/peptides_demo.csv')
    ap.add_argument('--features_root', default='data/features/replicates')
    ap.add_argument('--out_dir', default='data/models')
    args = ap.parse_args()

    y_map = load_labels(args.labels_csv)
    X = load_feature_set(args.features_root)
    X = X[X['id'].isin(y_map.keys())]
    y = X['id'].map(y_map)
    X_feat = X[['insertion_nm','contacts_6A']]

    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('xgb', XGBRegressor(n_estimators=200, max_depth=3, learning_rate=0.1, subsample=0.9, colsample_bytree=0.9, random_state=42))
    ])

    Xtr, Xte, ytr, yte = train_test_split(X_feat, y, test_size=0.33, random_state=42)
    pipeline.fit(Xtr, ytr)
    score = pipeline.score(Xte, yte)

    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, 'metrics.txt'), 'w') as f:
        f.write(f"R2: {score:.3f}\n")
    print(f"[OK] ML train abgeschlossen. R2={score:.3f}")


if __name__ == '__main__':
    main()
