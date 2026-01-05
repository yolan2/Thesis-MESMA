import pandas as pd

df = pd.read_csv('mesma_rf_comparison.csv')
total = len(df)
mesma_correct = (df['mesma_frac'] > 0.5).sum()
rf_correct = (df['rf_frac'] > 0.5).sum()
mesma_pct = mesma_correct / total * 100
rf_pct = rf_correct / total * 100
print(f'MESMA correctly predicted: {mesma_pct:.2f}% ({mesma_correct}/{total})')
print(f'RF correctly predicted: {rf_pct:.2f}% ({rf_correct}/{total})')
print('--- OVERALL FIT ---')
print(f'Mean Absolute Difference (MESMA vs RF): {df["diff_abs"].mean():.4f}')
print(f'RMSE (MESMA vs RF): {(df["diff_abs"] ** 2).mean() ** 0.5:.4f}')
print(f'Number of comparisons: {total}')