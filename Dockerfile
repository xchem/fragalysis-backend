FROM informaticsmatters/rdkit-python3-debian:latest
ENV PYTHONUNBUFFERED 1

USER root

WORKDIR /code
COPY requirements.txt /code/
RUN pip install -r requirements.txt

ADD . /code/
RUN apt-get --allow-releaseinfo-change update -y
RUN apt-get install -y nginx curl git default-libmysqlclient-dev redis-server
RUN apt-get install -y pandoc texlive-latex-base texlive-fonts-recommended
RUN chmod +x launch-stack.sh
COPY django_nginx.conf /etc/nginx/sites-available/default.conf
RUN ln -s /etc/nginx/sites-available/default.conf /etc/nginx/sites-enabled
COPY nginx.conf /etc/nginx/nginx.conf
RUN mkdir /srv/logs/
