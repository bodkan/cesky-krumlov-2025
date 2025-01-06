all: slides handouts exercises

slides:
	sed -i '/### handouts/d' slides.qmd
	sed -i '/### slides/ s/^# //' slides.qmd
	quarto publish quarto-pub --id b9e6d690-356c-4137-a475-4518fea71bbb slides.qmd
	git checkout slides.qmd
	git add slides*; git commit -m "Update slides"; git push

handouts:
	sed -i '/### slides/d' slides.qmd
	sed -i '/### handouts/ s/^# //' slides.qmd
	quarto publish quarto-pub --id cf112aa3-787f-4d8d-89b3-06b256aa0fbd slides.qmd
	git checkout slides.qmd
	git add slides*; git commit -m "Update handouts"; git push

exercises:
	quarto publish quarto-pub --id 21452ae4-076d-46a3-a620-c85e064edb1e exercises.qmd
	git add exercises*; git commit -m "Update exercises"; git push
