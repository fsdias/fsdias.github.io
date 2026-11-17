# To run in a terminal (not the R console), from this folder.
# The site is published to GitHub Pages at https://fsdias.github.io/

# 1. Save the source to GitHub (branch main)
git add -A
git commit -m "Update website"
git push

# 2. Render the site and publish it (branch gh-pages)
quarto publish gh-pages --no-prompt
