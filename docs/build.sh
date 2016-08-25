#!/usr/bin/env bash

rm -rf site/*

mkdocs build

rm -rf __pycache__
rm -rf featuredArticles
rm -rf css
rm -rf documentation
rm -rf fonts
rm -rf img
rm -rf js
rm -rf mkdocs

mv site/* ./

rm -rf site

rm -rf __pycache__
rm -rf __init__*
