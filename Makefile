all: slides handouts exercises

slides:
	sed -i '/### handouts/ s/^/#/' slides.qmd
	quarto render slides.qmd -o slides.html
	quarto publish quarto-pub --id 8675b273-c644-465d-98d3-08bf4963beb6 slides.qmd
	sed -i '/### handouts/ s/^#//' slides.qmd
	git add slides.qmd; git commit -m "Update slides.qmd"; git push
	git add slides.html; git commit -m "Update slides.html"; git push

handouts:
	sed -i '/### slides/ s/^/#/' slides.qmd
	quarto render slides.qmd -o handouts.html
	quarto publish quarto-pub --id 5eb6f220-2b2b-4843-b6f6-bbfe8fd3fe6f slides.qmd
	sed -i '/### slides/ s/^#//' slides.qmd
	git checkout slides.html
	git add slides.qmd; git commit -m "Update slides.qmd"; git push
	git add handouts.html; git commit -m "Update handouts.html"; git push

exercises:
	quarto publish quarto-pub --id 31de4a78-15a1-4cb8-a3ac-296a0bbc7d8b exercises.qmd
	git add exercises*; git commit -m "Update exercises"; git push
