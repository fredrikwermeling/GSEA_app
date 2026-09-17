FROM docker.io/library/httpd:alpine

COPY --chown=www-data:www-data . /usr/local/apache2/htdocs
