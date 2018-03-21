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
sendPage="sendPage_dir.php"

update:
	@make -s -j send_all_page
	@php $(sendPage) $(pages_dir)/navigation $(api_url) $(U) $(P)

update_old:
	@make -s send_genes
	@make -s send_metabolites
	@make -s send_pathways
	@make -s send_reactions
	@php $(sendPage) $(pages_dir)/navigation $(api_url) $(U) $(P)

send_all_page: send_genes send_metabolites send_pathways send_reactions


send_genes:
	@php $(sendPage) $(pages_dir)/genes $(api_url) $(U) $(P)

send_metabolites:
	@php $(sendPage) $(pages_dir)/metabolites $(api_url) $(U) $(P)

send_pathways:
	@php $(sendPage) $(pages_dir)/pathways $(api_url) $(U) $(P)

send_reactions:
	@php $(sendPage) $(pages_dir)/reactions $(api_url) $(U) $(P)



