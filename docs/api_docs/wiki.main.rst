===================
Wiki
===================
Description:
Wiki contains all the required ressources to create wiki (based on `mediaWiki <https://www.mediawiki.org/wiki/MediaWiki/>`_) for metabolic network.
A wiki can be created within a docker container or can be deployed on a server.


Docker version
==============
The docker version is by default made to be used on local machine, the wiki create are not protected, with default admin account, thus for update automatization.

It is possible to use this version for remote machine, but is highly recommended to change the password of the sql database and the default admin account for wiki.

Requirement
-----------
	a. php5-fpm

	b. docker

Get docker image
----------------

The image can be obtained with the command::

	sudo docker pull dyliss/wiki-img

or from sources::
	
	cd padmet-utils/wiki/docker_image_data
	sudo docker build -t dyliss/wiki-img .

Run docker container
--------------------

::

	sudo docker run -d -p 8888:80 -v /absolut/path/to/your/workdir:/shared --name 'wiki' dyliss/wiki-img

-p: set port, if 80 already used, change it.

-v: used to mount a folder in the container. Used to share folder of wiki-pages required for the wiki update

Wiki deployement inside container
---------------------------------

.. note:: It is possible to run the next commands from the container directly. To do this, it is required to go inside the container with the command:
	::

		sudo docker exec -ti wiki bash
	If done, all the next commands will be run without 'sudo docker exec wiki'

1. **Create a new wiki**::

	sudo docker exec wiki wiki --init='wiki_name'
	wiki --init='wiki_name' (if inside container)
	

Follow the instruction on the terminal, this should be prompted::

	Checking wiki id 'wiki_name'...
	Checking wiki id 'wiki_name': OK
	Checking if the prefix 'wiki_name_' is already used in the database...
	Checking the if the prefix 'wiki_name_' is already used in the database: OK
	Wiki initialization...
		Copying wiki folder
		Updating var in LocalSettings.php


	##############################################################
	MANUAL SETUP IS NOW REQUIRED. Access to this link from your browser:
		http://localhost:8888/wiki_name/mw-config/index.php
	Follow this instructions to setup mediawiki:
	Language:
		Continue->
	Existing wiki:
		Upgrade key: e94b13c7ca0b922f
		Continue->
	Welcome to MediaWiki!:
		Continue->
	Database settings:
		Continue->
	Name:
		Name of wiki: metabolic_network
		Administrator account: /!\ Use exactly the same to allow the bot to upload the pages automatically
			Your username: admin
			Password: 123456789
		I'm bored already, just install the wiki.
		Continue->
	Install:
		Continue->
		Continue->
		Do not save the LocalSettings file
	##############################################################
	When the previous setup is done, press enter to continue...

.. note:: If the name of the wiki is already used, it is mandatory to choose an other or remove the wiki using this name.

Wiki access
-----------
If the setup is complet, the wiki should be available on the link: 'http://localhost:8888/wiki_name/index.php'

wiki folder inside the container is in: /var/www/html/

Upload the wikipage of your metabolic network
---------------------------------------------

From your metabolic network you need to create wiki-pages with `padmet-utils <https://github.com/AuReMe/padmet-utils>`_ script, see `Documentation <https://padmet-utils.readthedocs.io/en/latest/api_docs/scripts.connection.html#module-scripts.connection.wikiGenerator>`_

Once you have your folder of wiki-pages, put this folder in your workdir specified during the creation of the container, **'/absolut/path/to/your/workdir'**.

The final step is to update the wiki from the folder of wiki-pages::

	sudo docker exec wiki wiki --wiki='wiki_name' --update=/shared/folder_name

	wiki --wiki='wiki_name' --update=/shared/folder_name (from the container)

During the wiki setup, if the administrator account is not the default one (admin, 123456789), it is required to add the args --username='XXX' --password='XXX'

If the wiki to update is on an remote machine, it is required to add arg --api_url='your/wiki/url/api.php'

Clean a wiki, delete all pages
------------------------------
The wiki will still be available but without any previous pages.

::

	sudo docker exec wiki wiki --wiki='wiki_name' --clean

	wiki --wiki='wiki_name' --clean (from the container)

Remove a wiki and clean sql database associated
-----------------------------------------------

::

	sudo docker exec wiki wiki --wiki='wiki_name' --remove

	wiki --wiki='wiki_name' --remove (from the container)


