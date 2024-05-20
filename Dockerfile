FROM python:3.11.9-slim-bullseye  AS python-base

ENV PYTHONUNBUFFERED 1
ENV PYTHONDONTWRITEBYTECODE 1

USER root

RUN apt-get update -y && \
    apt-get install --no-install-recommends -y \
      default-libmysqlclient-dev \
      nginx \
      pandoc \
      texlive-latex-base \
      texlive-latex-recommended \
      lmodern \
      texlive-fonts-recommended && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# another stage for poetry installation. this ensures poetry won't end
# up in final image where it's not needed
FROM python-base AS poetry-base

ARG POETRY_VERSION=1.7.1
RUN pip install --no-cache-dir poetry==${POETRY_VERSION}

WORKDIR /
COPY poetry.lock pyproject.toml /

# POETRY_VIRTUALENVS_IN_PROJECT tells poetry to create the venv to
# project's directory (.venv). This way the location is predictable
RUN POETRY_VIRTUALENVS_IN_PROJECT=true poetry install --no-root --only main --no-directory

# final stage. only copy the venv with installed packages and point
# paths to it
FROM python-base as final

COPY --from=poetry-base /.venv /.venv

ENV PYTHONPATH="${PYTHONPATH}:/.venv/lib/python3.11/site-packages/"
ENV PATH=/.venv/bin:$PATH

WORKDIR /srv/logs
WORKDIR /code/logs
WORKDIR /code

COPY nginx.conf /etc/nginx/nginx.conf
COPY django_nginx.conf /etc/nginx/sites-available/default.conf
COPY proxy_params /etc/nginx/frag_proxy_params
RUN ln -s /etc/nginx/sites-available/default.conf /etc/nginx/sites-enabled

COPY . ./
