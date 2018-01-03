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
MAINTAINER "Meziane A."
LABEL Version="0.4"
LABEL Vendor="INRIA-Dyliss team"
LABEL Description="MediaWiki installation in Docker container that runs on Nginx, MySQL, php5-fpm"

# Settings
ENV MYSQL_DB_ROOT_PWD="db_wiki_pwd"
ENV MEDIA_WIKI_ARCHIVE="mediawiki-1.26.2.tar.gz"

# Update repositories & install packages
RUN export DEBIAN_FRONTEND="noninteractive"
RUN apt-get update;\
    apt-get install -y wget nginx php5-fpm php5-mysql php5-apcu expect supervisor apt-utils;\
    echo mysql-server mysql-server/root_password password $MYSQL_DB_ROOT_PWD | debconf-set-selections;\
    echo mysql-server mysql-server/root_password_again password $MYSQL_DB_ROOT_PWD | debconf-set-selections;\
    apt-get install -y mysql-server;\
    apt-get -y -q clean;\
    apt-get -y -q autoremove;\
    apt-get install -y -q curl;\
    apt-get install -y -q vim;\
    apt-get install -y -q php5-curl;\
    apt-get install -y -q python-pip;\
    apt-get install -y -q imagemagick;\
    pip install docopt

# MediaWiki install
COPY wiki_template /home/wiki_template
COPY forms /home/forms
COPY wikiManager /home/wikiManager
COPY wiki /bin/
RUN chmod +x /bin/wiki
RUN mkdir /shared
RUN echo 'alias DB="mysql --user=root --password=$MYSQL_DB_ROOT_PWD"' >> ~/.bashrc
# Config Nginx
COPY default /etc/nginx/sites-available/

# Config php5-fpm upload size
RUN sed -i -e 's#upload_max_filesize = 2M#upload_max_filesize = 25M#g' /etc/php5/fpm/php.ini
RUN sed -i -e 's#post_max_size = 8M#post_max_size = 25M#g' /etc/php5/fpm/php.ini

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
