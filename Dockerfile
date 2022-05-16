FROM informaticsmatters/rdkit-python3-debian:Release_2021_09_2
ENV PYTHONUNBUFFERED 1

USER root

WORKDIR /code
COPY requirements.txt /code/
RUN pip install --upgrade pip && \
    pip install -r requirements.txt

ADD . /code/
RUN apt-get --allow-releaseinfo-change update -y && \
    apt-get install -y \
      curl \
      default-libmysqlclient-dev \
      git \
      nginx \
      pandoc \
      redis-server \
      texlive-latex-base \
      texlive-fonts-recommended && \
    chmod +x launch-stack.sh

COPY django_nginx.conf /etc/nginx/sites-available/default.conf
RUN ln -s /etc/nginx/sites-available/default.conf /etc/nginx/sites-enabled
COPY nginx.conf /etc/nginx/nginx.conf
WORKDIR /srv/logs/

WORKDIR /code
