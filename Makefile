all: slides handouts exercises

slides:
	sed -i '/### handouts/d' slides.qmd
	quarto render -o slides.html slides.qmd
	quarto publish quarto-pub --id 27532d2d-db64-47ed-a938-9a4102c43a6a slides.qmd
	git checkout slides.qmd
	git add slides*; git commit -m "Update slides"; git push

handouts:
	sed -i '/### slides/d' slides.qmd
	sed -i '/### handouts/ s/^# //' slides.qmd
	quarto render -o handouts.html slides.qmd
	quarto publish quarto-pub --id e1b4dcad-9c7b-40ca-bfc0-47cfbdc91bed slides.qmd
	git checkout slides.qmd
	git add slides*; git commit -m "Update handouts"; git push

exercises:
	quarto render -o exercises.html exercises.qmd
	quarto publish quarto-pub --id 21452ae4-076d-46a3-a620-c85e064edb1e exercises.qmd
	git add exercises*; git commit -m "Update exercises"; git push
