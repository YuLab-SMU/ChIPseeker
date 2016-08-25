build: feature index
	sh build.sh

feature:
	cd scripts;\
	Rscript featured_article.R

#index:
#	cd scripts;\
#	Rscript make_index.R

index:
	cd docs;\
	Rscript -e 'rmarkdown::render("../private/index.rmd", "md_document")';\
	mv ../private/index.md ./
