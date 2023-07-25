FROM python:3.10.12

ENV PYTHONUNBUFFERED 1
ENV PYTHONDONTWRITEBYTECODE 1
ENV POETRY_VERSION=1.5.1
ENV POETRY_HOME=/opt/poetry
ENV PATH="${PATH}:${POETRY_HOME}/bin"

USER root

# Install required packages.
# gnupg (and bzip2 & wget) is actually used by the stack,
# we load them here to simplify the stack's Dockerfile.
RUN apt-get update --allow-insecure-repositories -y && \
    apt-get install --no-install-recommends -y \
      default-libmysqlclient-dev \
      gnupg \
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

COPY poetry.lock pyproject.toml ./

RUN curl -sSL https://install.python-poetry.org | python - && \
    poetry export -f requirements.txt --without dev --output requirements.txt && \
    pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir --requirement requirements.txt

COPY . ./
