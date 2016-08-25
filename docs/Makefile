push: build
	git add .; git commit -m 'docs'; git push -u origin gh-pages

build: feature
	sh build.sh

feature:
	cd scripts;\
	Rscript featured_article.R

