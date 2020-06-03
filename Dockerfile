FROM informaticsmatters/rdkit-python-debian:latest
ENV PYTHONUNBUFFERED 1
ADD . /code/
WORKDIR /code
USER root
RUN apt-get --allow-releaseinfo-change update -y
RUN apt-get install -y nginx curl git default-libmysqlclient-dev redis-server
RUN pip install -r requirements.txt
RUN chmod +x launch-stack.sh
COPY django_nginx.conf /etc/nginx/sites-available/default.conf
RUN ln -s /etc/nginx/sites-available/default.conf /etc/nginx/sites-enabled
COPY nginx.conf /etc/nginx/nginx.conf
RUN mkdir /srv/logs/