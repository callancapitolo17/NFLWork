# EfficiencyLandscape

Interactive Shiny app for exploring NFL offensive efficiency landscapes. Shows team efficiency as scatter plots with EPA and success rate metrics.

## Files

| File | Description |
|------|-------------|
| `app.R` | Shiny app with season selector (2000-2025). Loads play-by-play data, calculates offensive efficiency, renders plotly scatter plot with team logos. |
| `49ersOffense.png` | Pre-rendered example visualization |

## Running Locally

```r
shiny::runApp("EfficiencyLandscape")
```

## Dependencies

- `shiny`, `nflfastR`, `plotly`, `tidyverse`
