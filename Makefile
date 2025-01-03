slides:
	rm -r slides_files
	quarto publish quarto-pub --id 98b8709c-8190-4cdc-900b-4df3953ec283 slides.qmd
	git add slides*

onepage:
	rm -r slides_files
	quarto publish quarto-pub --id 7f5d556f-ca21-4c53-b06c-693befbfa994 slides.qmd
