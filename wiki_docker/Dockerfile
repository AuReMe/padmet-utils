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



# Dockerfile is located in the root of the context.
# Warning: Do not use your root directory, /, 
# as the PATH as it causes the build to transfer the entire contents of your hard drive to the Docker daemon.


FROM debian:jessie
MAINTAINER "Pierre V. - Meziane A."
LABEL Version="0.3"
LABEL Vendor="INRIA-Dyliss team"
LABEL Description="MediaWiki installation in Docker container that runs on Nginx, MySQL, php5-fpm"

# Settings
ENV MYSQL_DB_ROOT_PWD="db_wiki_pwd"
ENV MEDIA_WIKI_ARCHIVE="mediawiki-1.26.2.tar.gz"
ENV NGINX_DIR="/var/www/html/wiki"

# Update repositories & install packages
RUN export DEBIAN_FRONTEND="noninteractive"
RUN apt-get update;\
    apt-get install -y wget nginx php5-fpm php5-mysql php5-apcu expect supervisor apt-utils;\
    echo mysql-server mysql-server/root_password password $MYSQL_DB_ROOT_PWD | debconf-set-selections;\
    echo mysql-server mysql-server/root_password_again password $MYSQL_DB_ROOT_PWD | debconf-set-selections;\
    apt-get install -y mysql-server;\
    apt-get -y -q clean;\
    apt-get -y -q autoremove;\
    apt-get install -y curl;\
    apt-get install -y vim

# MediaWiki install
# Export MediaWiki into Nginx dir & auto strip 'mediawiki-1.26.2' root directory
RUN cd /root/;\
	mkdir $NGINX_DIR;\
	wget https://releases.wikimedia.org/mediawiki/1.26/$MEDIA_WIKI_ARCHIVE;\
	tar -xzf $MEDIA_WIKI_ARCHIVE -C $NGINX_DIR --strip-components=1;\
	rm $MEDIA_WIKI_ARCHIVE

# Maybe some bugs in AUFS => switch from aufs to devicemapper or btrfs to resolve them.
#RUN mv /root/mediawiki-* $NGINX_DIR


# Config MediaWiki
COPY LocalSettings.php $NGINX_DIR/LocalSettings.php
RUN rm /var/www/html/index.nginx-debian.html;\
	cd $NGINX_DIR ;\
	curl -sS https://getcomposer.org/installer | php;\
	php composer.phar require mediawiki/semantic-media-wiki "2.4.5";\
	php composer.phar update

# Config Nginx
COPY default /etc/nginx/sites-available/

# Config php5-fpm upload size
RUN sed -i -e 's#upload_max_filesize = 2M#upload_max_filesize = 25M#g' /etc/php5/fpm/php.ini

# Change MySQL root password (don't forget rupture in simple quotes for the interpretation of variable)
RUN sed -i -e 's#$wgDBpassword = "";#$wgDBpassword = "'$MYSQL_DB_ROOT_PWD'";#g' $NGINX_DIR/LocalSettings.php

# Config MySQL
RUN service mysql start
RUN expect << END_OF_EXPECT ;\
        set timeout 1\
\
        spawn mysql_secure_installation\
\
        expect {"Enter current password for root (enter for none):"}\
        send "$MYSQL_DB_ROOT_PWD\n"\
\
        expect {Change the root password? [Y/n]}\
        send "n\r"\
\
        expect {Remove anonymous users? [Y/n]}\
        send "y\r"\
\
        expect {Disallow root login remotely? [Y/n]}\
        send "y\r"\
\
        expect {Remove test database and access to it? [Y/n]}\
        send "y\r"\
\
        expect {Reload privilege tables now? [Y/n]}\
        send "y\r"\
\
        expect eof\
END_OF_EXPECT\

# Restore SQL dump
# ADD tar.gz /root/
# mysql --user=root --password=$MYSQL_DB_ROOT_PWD my_wiki < tut_backup.sql

# Informs Docker that the container listens on the specified network ports at runtime. 
# EXPOSE does not make the ports of the container accessible to the host
EXPOSE 80 3306

# Command for an 'executing container'
# Think to Python supervisord because the following line is disgusting
#CMD service nginx start & service mysql start & service php5-fpm start && tail -F /var/log/nginx/error.log
COPY supervisord.conf /etc/supervisor/conf.d/
CMD ["/usr/bin/supervisord", "-c", "/etc/supervisor/supervisord.conf"]
