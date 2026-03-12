import sys

with open("fit_veg_mixture_mesma.R", "r", encoding="utf-8") as f:
    text = f.read()

old_y_target = """      y_target <- y
      y_target[!is.finite(y_target)] <- 0"""

new_y_target = """      y_target <- y
      # y_target[!is.finite(y_target)] <- 0"""

text = text.replace(old_y_target, new_y_target)

with open("fit_veg_mixture_mesma.R", "w", encoding="utf-8") as f:
    f.write(text)

print("Patch 2 applied")
