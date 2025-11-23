# Thesis MESMA

This repository contains the code and data flow for the Master's thesis project implementing a multiple endmember spectral mixture analysis (MESMA) workflow for vegetation phenology mapping.

Structure:
- fit_veg_mixture_mesma.R — main MESMA fit script
- transform_phenology.py — transformer used for phenology data
- Other helper scripts and resources

I created a local git repository and added a .gitignore including modis, landsat, phenology_results and virtualenv patterns.

To push this repository to GitHub (example):

1. Create a new repo named `Thesis MESMA` on GitHub.
2. Then run (PowerShell):

```powershell
git remote add origin https://github.com/<your-username>/Thesis-MESMA.git
git branch -M main
git push -u origin main
```

If you'd like me to create the remote GitHub repository now, tell me your GitHub preferences and I can help (but I can’t push without your credentials / an access token).
