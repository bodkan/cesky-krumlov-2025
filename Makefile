all: slides handouts exercises

slides:
	sed -i '/### handouts/ s/^/#/' slides.qmd
	quarto render slides.qmd -o slides.html
	quarto publish quarto-pub --id 27532d2d-db64-47ed-a938-9a4102c43a6a slides.qmd
	sed -i '/### handouts/ s/^#//' slides.qmd
	git add slides*; git commit -m "Update slides"; git push

handouts:
	sed -i '/### slides/ s/^/#/' slides.qmd
	quarto render slides.qmd -o handouts.html
	quarto publish quarto-pub --id e1b4dcad-9c7b-40ca-bfc0-47cfbdc91bed slides.qmd
	sed -i '/### slides/ s/^#//' slides.qmd
	git add handouts.html; git commit -m "Update handouts"; git push

exercises:
	quarto render exercises.qmd -o exercises.html
	quarto publish quarto-pub --id 21452ae4-076d-46a3-a620-c85e064edb1e exercises.qmd
	git add exercises*; git commit -m "Update exercises"; git push
