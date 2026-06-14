# Deploying HTML reports with GitHub Pages

This repository uses the `docs/` folder as the GitHub Pages site root.

## 1. Confirm the expected structure

```text
docs/
├── index.html
├── DEPLOYMENT.md
├── reports/
│   ├── integrated_tb_amr_mtbc_phylogenomics_report.html
│   └── tbprofiler_combined_report.html
└── assets/
    └── .gitkeep
```

## 2. Copy reports into `docs/reports/`

After a successful pipeline run:

```bash
mkdir -p docs/reports docs/assets

cp /path/to/integrated_tb_amr_mtbc_phylogenomics_report.html docs/reports/
cp /path/to/tbprofiler_combined_report.html docs/reports/
```

If the reports depend on images, CSS, JavaScript, or other assets, copy those into `docs/assets/` and ensure the HTML uses relative paths.

```bash
cp -R /path/to/assets/* docs/assets/
```

## 3. Commit and push

```bash
git add docs/
git commit -m "Update GitHub Pages reports"
git push origin main
```

## 4. Enable GitHub Pages

In the GitHub repository:

1. Open **Settings**.
2. Select **Pages**.
3. Under **Build and deployment**, choose **Deploy from a branch**.
4. Select:

```text
Branch: main
Folder: /docs
```

5. Save.

## 5. Open the deployed site

```text
https://gmboowa.github.io/TB-AMR-MTBC-Phylogenomics/
```

## 6. Troubleshooting

| Problem | Likely cause | Fix |
|---|---|---|
| Page gives 404 | GitHub Pages not enabled or still building | Wait a few minutes and confirm Pages settings |
| Cards open broken links | Report filenames differ from `index.html` links | Rename reports or update `docs/index.html` |
| Report opens but images are missing | Assets were not copied | Copy assets into `docs/assets/` and use relative paths |
| GitHub Pages uses old content | Browser cache or old commit | Hard refresh and confirm latest commit is on `main` |
