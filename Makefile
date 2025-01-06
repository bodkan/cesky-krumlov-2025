all: slides handouts exercises

slides:
	sed -i '/### handouts/d' slides.qmd
	sed -i '/### slides/ s/^  #//' slides.qmd
	quarto publish quarto-pub --id b9e6d690-356c-4137-a475-4518fea71bbb slides.qmd
	git checkout slides.qmd
	git add slides*; git commit -m "Update slides"; git push

handouts:
	sed -i '/### slides/d' slides.qmd
	sed -i '/### handouts/ s/^  #//' slides.qmd
	quarto publish quarto-pub --id 7f5d556f-ca21-4c53-b06c-693befbfa994 slides.qmd
	git checkout slides.qmd
	git add slides*; git commit -m "Update handouts"; git push

exercises:
	quarto publish quarto-pub --id 7f5d556f-ca21-4c53-b06c-693befbfa994 slides.qmd
	git add exercises*; git commit -m "Update exercises"; git push
