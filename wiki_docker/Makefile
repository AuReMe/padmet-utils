#This file is part of padmet-utils.
#
#padmet-utils is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#padmet-utils is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with padmet-utils. If not, see <http://www.gnu.org/licenses/>.

#Makefile for media_wiki Docker container
name=default_wiki
sendPage="wikiManager/sendPage_dir.php"

all: build run

b: build

r: run

h: help

l: list

i: init

build:
	docker build -t dyliss/media_wiki_docker .

run:
	@docker run -d -p 8888:80 --name $(name) dyliss/media_wiki_docker
	@echo "Follow the README to configure the wiki"
	@echo "http://localhost:8888/wiki/mw-config/index.php"
	@echo "upgrade key: 0880a534700be3f7"
	@echo "Default admin account: Username=admin; Pwd=123456789"

init:
	@docker exec $(name) /var/www/html/wiki/maintenance/update.php
	@make -s -j send_all_page
	@php $(sendPage) $(pages_dir)/navigation
	@echo "Setup done. Use this URL to access to the wiki:"
	@echo "http://localhost:8888/wiki/index.php?title=Main_Page"

send_all_page: send_genes send_metabolites send_pathways send_reactions


send_genes:
	@php $(sendPage) $(pages_dir)/genes

send_metabolites:
	@php $(sendPage) $(pages_dir)/metabolites

send_pathways:
	@php $(sendPage) $(pages_dir)/pathways

send_reactions:
	@php $(sendPage) $(pages_dir)/reactions


help:
	@echo "Command:	Result:"
	@echo "<make build>	Build the image."
	@echo "<make run>	Run the image in new container, specify the name with the argument name=XXX. The default name is: "$(name)"."
	@echo "<make init>	Initialization of mediawiki after finishing the setup, specify the name with the argument name=XXX. The default name is: "$(name)". Also specify the directory with all the wiki pages to send with the argument pages_dir=XXX. No default directory."
	@echo "<make help>	Display this help."
	@echo "<make version>	Display informations about the image."
	@echo "Once the image is built and the container is up, you can access the web server through your browser at http://127.0.0.1:8888"

version:
	@docker inspect --format='Description: 	{{.Config.Labels.Description}}' inria/media_wiki_docker
	@docker inspect --format='Vendor: 	{{.Config.Labels.Vendor}}' inria/media_wiki_docker
	@docker inspect --format='Authors: 	{{.Author}}' inria/media_wiki_docker
	@docker inspect --format='Version: 	{{.Config.Labels.Version}}' inria/media_wiki_docker

	@docker inspect --format='DockerVersion:	{{.DockerVersion}}' inria/media_wiki_docker
	@docker inspect --format='Architecture:	{{.Architecture}}' inria/media_wiki_docker
	@docker inspect --format='OS: 		{{.Os}}' inria/media_wiki_docker
	@docker inspect --format='Size: 		{{.Size}} bytes' inria/media_wiki_docker


