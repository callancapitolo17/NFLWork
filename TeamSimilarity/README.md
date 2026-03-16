# TeamSimilarity

Interactive Shiny app that finds historically similar NFL teams using offensive and defensive scouting profiles with DBSCAN clustering. Deployed on shinyapps.io.

## Files

| File | Description |
|------|-------------|
| `app.R` | Shiny app that loads 2006-2024 play-by-play data, builds comprehensive offensive/defensive profiles, applies DBSCAN clustering to identify similar teams across eras, and displays results in formatted gt() tables. |
| `rsconnect/` | Deployment configuration for shinyapps.io/callan-capitolo |

## Running Locally

```r
shiny::runApp("TeamSimilarity")
```

## Dependencies

- `shiny`, `nflfastR`, `dbscan`, `gt`, `gtExtras`, `tidyverse`
