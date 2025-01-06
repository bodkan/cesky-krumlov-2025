all: slides onepage exercises

slides:
	cp slides.qmd slides.qmd.bak
	sed -i '/### onepage/d' slides.qmd
	quarto publish quarto-pub --id b9e6d690-356c-4137-a475-4518fea71bbb slides.qmd
	mv slides.qmd.bak slides.qmd
	git add slides*; git commit -m "Update slides"; git push

onepage:
	cp slides.qmd slides.qmd.bak
	sed -i '/### slides/d' slides.qmd
	sed -i '/### onepage/ s/^  #//' slides.qmd
	quarto publish quarto-pub --id 7f5d556f-ca21-4c53-b06c-693befbfa994 slides.qmd
	mv slides.qmd.bak slides.qmd
	git add slides*; git commit -m "Update handouts"; git push

exercises:
	quarto publish quarto-pub --id 7f5d556f-ca21-4c53-b06c-693befbfa994 slides.qmd
	git add exercises*; git commit -m "Update exercises"; git push
