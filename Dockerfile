FROM informaticsmatters/rdkit-python3-debian:Release_2021_09_2
ENV PYTHONUNBUFFERED 1

USER root

# Install required packages.
# bzip2, gnupg & wget are actually used by the stack,
# we load them here to simplify the stack Dockerfile.
RUN apt-get --allow-releaseinfo-change update -y && \
    apt-get install -y \
      bzip2 \
      curl \
      default-libmysqlclient-dev \
      git \
      gnupg \
      nginx \
      pandoc \
      redis-server \
      texlive-latex-base \
      texlive-fonts-recommended \
      wget

COPY django_nginx.conf /etc/nginx/sites-available/default.conf
RUN ln -s /etc/nginx/sites-available/default.conf /etc/nginx/sites-enabled
COPY nginx.conf /etc/nginx/nginx.conf

WORKDIR /code
ADD . /code/
RUN chmod 755 launch-stack.sh && \
    chmod 755 makemigrations.sh && \
    chmod 755 docker-entrypoint.sh

RUN pip install --upgrade pip && \
    pip install -r requirements.txt

WORKDIR /srv/logs
WORKDIR /code/logs
WORKDIR /code
