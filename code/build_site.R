# Load workflowr package
library(workflowr)

# Build the website
wflow_build()

# View the website
wflow_view()

# Publish the website (commits and pushes to GitHub)
# wflow_publish(c("analysis/index.Rmd", "analysis/about.Rmd", "analysis/license.Rmd"),
#               "Initial workflowr setup")

# Check workflowr status
wflow_status()
