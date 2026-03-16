# EPAbyWeek

Interactive Shiny app for tracking cumulative EPA (Expected Points Added) across NFL weeks. Deployed on shinyapps.io.

## Files

| File | Description |
|------|-------------|
| `EPAbyWeek.R` | Shiny app with selectInputs for X/Y axes (Offense/Defense, yards_gained/success/epa). Renders plotly scatter plots showing cumulative team performance as the season progresses. |
| `rsconnect/` | Deployment configuration for shinyapps.io/callan-capitolo |

## Running Locally

```r
shiny::runApp("EPAbyWeek")
```

## Dependencies

- `shiny`, `nflfastR`, `plotly`, `tidyverse`
