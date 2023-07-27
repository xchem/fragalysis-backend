FROM python:3.11.4-slim-bullseye

# The build process relies on the presence of a requirements.txt file in project root.
# This is achived prior to a container build using poetry
# and is a responsibility of the developer or CI process: -
#
#   poetry export --without-hashes --without dev --output requirements.txt

ENV PYTHONUNBUFFERED 1
ENV PYTHONDONTWRITEBYTECODE 1

USER root

RUN apt-get update --allow-insecure-repositories -y && \
    apt-get install --no-install-recommends -y \
      default-libmysqlclient-dev \
      nginx \
      pandoc \
      texlive-latex-base \
      texlive-fonts-recommended && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY django_nginx.conf /etc/nginx/sites-available/default.conf
RUN ln -s /etc/nginx/sites-available/default.conf /etc/nginx/sites-enabled
COPY nginx.conf /etc/nginx/nginx.conf

WORKDIR /srv/logs
WORKDIR /code/logs
WORKDIR /code

COPY . ./

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir --requirement requirements.txt
