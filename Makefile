all: slides onepage

slides:
	cp slides.qmd slides.qmd.bak
	sed -i '/### onepage/d' slides.qmd
	quarto publish quarto-pub --id 98b8709c-8190-4cdc-900b-4df3953ec283 slides.qmd
	mv slides.qmd.bak slides.qmd
	git add slides*; git commit -m "Update slides"

onepage:
	cp slides.qmd slides.qmd.bak
	sed -i '/### slides/d' slides.qmd
	sed -i '/### onepage/ s/^  #//' slides.qmd
	quarto publish quarto-pub --id 7f5d556f-ca21-4c53-b06c-693befbfa994 slides.qmd
	mv slides.qmd.bak slides.qmd
	git add slides*; git commit -m "Update handouts"

exercises:
	quarto render exercises.qmd
	git add exercises*; git commit -m "Update exercises"
