FROM python:3.13.18-slim-trixie  AS python-base

ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

USER root

# Apply the Debian (trixie) security updates to the base image packages
# and upgrade pip (the base image's pip has known vulnerabilities).
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install --no-install-recommends -y \
      default-libmysqlclient-dev \
      nginx \
      pandoc \
      texlive-latex-base \
      texlive-latex-recommended \
      lmodern \
      pigz \
      p7zip-full \
      rsync \
      build-essential \
      texlive-fonts-recommended && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir --upgrade pip

# another stage for poetry installation. this ensures poetry won't end
# up in final image where it's not needed
FROM python-base AS poetry-base

ARG POETRY_VERSION=2.1.4
RUN pip install --no-cache-dir poetry==${POETRY_VERSION}

WORKDIR /
COPY poetry.lock pyproject.toml /

# POETRY_VIRTUALENVS_IN_PROJECT tells poetry to create the venv to
# project's directory (.venv). This way the location is predictable.
# Only the 'main' dependencies are installed - the unit tests run with pytest
# on the host / CI runner (see run-unit-tests.sh), not inside this image.
# Poetry does not manage the pip that virtualenv seeds into the venv,
# so upgrade it explicitly (the seeded pip has known vulnerabilities).
RUN POETRY_VIRTUALENVS_IN_PROJECT=true poetry install --no-root --only main --no-directory && \
    /.venv/bin/pip install --no-cache-dir --upgrade pip

# final stage. only copy the venv with installed packages and point
# paths to it
FROM python-base AS final

COPY --from=poetry-base /.venv /.venv

ENV PYTHONPATH="/.venv/lib/python3.13/site-packages/"
ENV PATH=/.venv/bin:$PATH
ENV TMPDIR=/tmp

WORKDIR /code/logs
WORKDIR /code

COPY nginx.conf /etc/nginx/nginx.conf
COPY django_nginx.conf /etc/nginx/sites-available/default.conf
COPY proxy_params /etc/nginx/frag_proxy_params
RUN ln -s /etc/nginx/sites-available/default.conf /etc/nginx/sites-enabled

COPY . ./

CMD ["./docker-entrypoint.sh"]
